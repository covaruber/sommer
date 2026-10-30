# Native GLMM epic: implementation record and remaining roadmap

Status reviewed: 2026-10-27, against the local working tree.

**Overall status: M1 correctness, M2 native modular/first-score core, and M3
agreed response-family set implemented; the bounded M4 inference, prediction,
and simulation implementation is in place; the first restricted M5
zero-inflated Poisson slice is implemented.** Broader scientific and package
qualification and the remaining M5 features are outstanding. Implemented does
not mean committed, release-ready, or validated for every covariance family
or platform.

### Milestone 1 execution update

Correctness changes implemented on 2026-10-26. Focused validation passed;
full release qualification was interrupted at the user's request to move to
M2. It is not a clean completed package check.

- Preparation forwards supplied contrasts and all missing-data policies.
  Prediction retains fitted coding, retained factor levels, transformed-term
  metadata, and a column-position map for matrix-valued terms.
- Observation masks are retained in the fit. Offset construction now directly
  uses the retained response length when no offset is supplied.
- Native Poisson validation rejects negative OR non-integral counts. Native
  precision symmetry, positive definiteness and determinant consistency are
  checked; numerical mode failures return explicit reasons.
- Outer evaluations are deterministic from a zero mode start; blanket
  `[-25, 25]` descriptor bounds have been removed. Descriptor transforms and
  covariance positive definiteness govern admissibility, with a
  machine-derived exponential-coordinate overflow/underflow guard.
- Fits distinguish optimizer convergence, conditional-mode convergence,
  curvature validity, boundary screening and availability of uncertainty.
  Invalid curvature and near-boundary fits warn and return unavailable
  standard errors, never clamped zero standard errors.
- Covariance evaluator/callback errors propagate. Expected numerical
  rejections return an infinite objective and record rejection reasons.
  Only the identified finite-difference Hessian failure is converted to
  unavailable uncertainty, with an explicit warning.
- Unused non-default Gaussian controls and unsupported prediction arguments
  are rejected. Working residuals are corrected; zero conditional variance
  or zero link derivative makes the corresponding residual undefined and
  produces an explicit error.
- Native scalar/multivariate likelihood and Hessian checks, repeated mode and
  covariance starts, filtering, contrasts/transformed terms, boundary,
  zero-parameter and failure regressions have been added.
- Neither Gaussian C++ engine body or objective has been changed.

**Correction to the initial gap assessment:** `obsInfo$included` is a logical
mask of the original row count, not a vector of retained integer indices.
Consequently the earlier no-offset allocation/subsetting was not itself a
demonstrated non-prefix-omission bug. It has been simplified and regression
tested, rather than claimed as a repaired observed failure.

## 1. Original target and architectural constraints

The work began with an audit of sommer's GLMM capabilities and a comparison
with ASReml, SAS GLIMMIX, and glmmTMB. The agreed direction was to:

- Preserve the existing PQL/IRLS route and Gaussian fitting behavior.
- Add a genuine approximate marginal GLMM likelihood, initially Laplace,
  permitting likelihood reporting and AIC/BIC.
- Retain sommer's structured random covariance and pedigree/genomic
  relationship-matrix strengths.
- Build a sommer-native, modular implementation without a fitting dependency
  on another mixed-model package.
- Make any disturbance to the Henderson and direct C++ engines explicit.
- Subsequently address distribution coverage, dispersion, zero inflation/
  hurdle models, prediction, inference, and diagnostics.

The initial implementation used an optional lme4 adapter. That approach was
rejected because it did not satisfy the native-engine requirement. It has been
replaced; it is not part of the current fitting architecture.

This document reconciles that target with the current source. It is not a
claim of feature parity with any competing package, nor a new exhaustive
literature review. The broad comparison identified several separate work
streams; completing a Laplace evaluator does not complete all of them.

## 2. What exists today

| Capability | Current status |
| --- | --- |
| PQL/IRLS and ordinary Gaussian routes | Retained; PQL remains the default GLMM estimator. |
| Native Laplace conditional likelihood | Binomial-logit/probit/cloglog, Poisson-log, restricted zero-inflated Poisson, NB2-log, Gaussian-identity and Gamma-log; grouped trials and mixed-family observation indexing connected. |
| Native random-effect conditional mode | Sparse Newton steps, backtracking, CHOLMOD factorization. |
| Laplace determinant correction | Implemented and covered by a scalar analytic regression test. |
| Outer optimization | R `stats::nlminb`, jointly over fixed effects and free covariance parameters; native full Laplace first scores, numerical uncertainty Hessian. |
| Shared model representation | Reuses sommer preparation, random design blocks, covariance descriptors, and relationship precision. |
| Random covariance coverage | Descriptor-native precision support is connected; tests cover scalar, unstructured, AR, random-slope, crossed and supplied relationship examples, including fixed parameters, not the entire covariance catalog. |
| Offsets | Formula offsets are connected; Poisson offset fitting is tested. |
| Likelihood output | `logLik`, AIC/BIC, fitted values, residuals, coefficients, summary/print. |
| Native shape/dispersion | NB2 shape, Gaussian/Gamma dispersion and the ZIP structural-zero probability estimated per family/trait by default; natural-scale output and valid-curvature delta-method SEs. |
| Prediction | Conditional, zero-random-effect, and restricted population-mean new-data predictions; post-fit `se.fit`/`cov.fit` are available within documented covariance and memory limits. |
| Post-fit uncertainty/PEV | Numerical-Hessian parameter uncertainty and sparse-panel or guarded full joint fixed/random PEV, computed only after fitting. |
| Main Gaussian C++ engine objectives | Unchanged. |
| Native GLMM covariance scores/updates | Full Laplace first scores connected to `nlminb`; Gaussian AI updates are not GLMM likelihood derivatives. Descriptor callbacks may themselves use numerical covariance derivatives. |
| Release qualification | Incomplete; see validation record and milestones below. |

### Public entry point

```r
fit_native <- mmes(
  outcome ~ treatment,
  random = ~ subject,
  data = observations,
  family = binomial(),
  control = mmes_control(
    task = list(
      objective = "ml",
      glmm = list(
        method = "laplace",
        maxit = 1000L,
        modeMaxit = 100L,
        modeTol = 1e-8
      )
    )
  )
)
```

`maxit` controls outer optimization. `modeMaxit` and `modeTol` control the
conditional random-effect mode solve. Gaussian AI iteration controls are not
substitutes for these controls.

### Current restrictions

Laplace requires Henderson formulation, ML, a full-rank fixed design and at
least one random term. Responses are one numeric vector (including
long-format mixed traits) or two-column binomial successes/failures.
Binomial successes/trials and count responses must be integers; Gamma
responses must be positive.

Not supported on the Laplace route:

- Families/links outside the declared M3 set and the restricted M5
  zero-inflated Poisson model, including beta-binomial, hurdle models and
  zero-inflated NB2.
- Structured residual covariance, supplied `W`, weights formulas, and `vcc`.
- Direct formulation and known-covariance solve mode. Laplace requires
  `output$inverse=0`; approximate fixed/random PEV is available after fitting
  through `postPEV()`.
- Gaussian engine accelerations, forced experimental representations, and
  other rejected non-default engine settings.
- General population integration, variance-component likelihood-ratio tests,
  and denominator-df Wald methods. Prediction standard errors are available
  for supported new-data conditional, zero-random-effect, and restricted
  population predictions; they remain plug-in curvature approximations.
- Zero-inflated Poisson is native-Laplace-only and has a separate scalar
  logit parameter by default, or a bounded fixed-effects zero formula through
  `zipoisson(zero=~covariates)`. Zero random effects remain unsupported.
  A zero-inflation probability can be confounded with the count intercept;
  broader identifiability and calibration qualification remains outstanding.

These restrictions do not describe the broader existing PQL feature set.
Rejecting a feature is not the same as implementing it.

## 3. Source map and actual engine disturbance

| Surface | Role in this epic |
| --- | --- |
| [R/FUN_glmm.R](../R/FUN_glmm.R) | `.mmes_laplace()` prepares the model, creates a native objective context, optimizes with native scores, and builds the fitted object. Also contains existing PQL. |
| [R/FUN_glmm_families.R](../R/FUN_glmm_families.R) | Native NB2 R family constructor, family specification/configuration, trial/response validation, observation indexing and fitted-family reporting helpers. |
| [R/mmes.R](../R/mmes.R) | Dispatches a Laplace request; private preparation-only contract exits after design/descriptor construction without Gaussian response preparation or fitting. |
| [mme_glmm_laplace.cpp](./mme_glmm_laplace.cpp) | Native likelihood, random-effect mode, fitting-local objective context and full first-score evaluation. |
| [mme_glmm_family.h](./mme_glmm_family.h) | Reusable stable family log density, mean, score, curvature and first/second curvature derivatives. |
| [mme_glmm_precision.h](./mme_glmm_precision.h) | Descriptor-native precision/determinant derivatives, cached relationship factorizations and sparse block assembly. |
| [mme_cholmod.h](./mme_cholmod.h) | Existing sparse factorization interface reused by the new module. |
| [R/FUN_mme_modules.R](../R/FUN_mme_modules.R) | GLMM controls, validation, and execution-plan reporting. |
| [R/FUN_mme_capabilities.R](../R/FUN_mme_capabilities.R) | Capability metadata. |
| [R/FUN_methods_mmes.R](../R/FUN_methods_mmes.R) | Likelihood, summary, print, fitted, residual, and inference guards. |
| [R/predict_mmes.R](../R/predict_mmes.R) | Laplace prediction branch. |
| [RcppExports.cpp](./RcppExports.cpp), [R/RcppExports.R](../R/RcppExports.R) | Generated binding for the native evaluator. |
| [tests/testthat/test-mmes-laplace.R](../tests/testthat/test-mmes-laplace.R) | Focused native regression tests. |
| [tests/testthat/test-mmes-laplace-families.R](../tests/testthat/test-mmes-laplace-families.R) | Independent expanded-family densities, derivatives, nuisance scores, grouped/mixed response and recovery tests. |
| [vignettes/sommer.glmm.Rmd](../vignettes/sommer.glmm.Rmd) | User-facing estimator description and limitations. |

**`ai_mme_sp2()` in [ai_mme_sp2.cpp](./ai_mme_sp2.cpp): no body or Gaussian
objective changes for the native Laplace implementation.**

**`ai_reml_direct_sp2()` in [mme_direct.cpp](./mme_direct.cpp): no body or
Gaussian objective changes; direct Laplace fitting remains unsupported.**

The preparation bridge calls `.mmes_fit()` with the actual response family
and private `.prepareOnly=TRUE` to obtain designs and descriptors. It exits
before Gaussian response preparation or fitting.

This is a shared **data/covariance/precision support seam**, coupled to a
separate native GLMM score and mode solver. It does not turn the Gaussian
AI/REML coefficient engines into GLMM optimizers.

## 4. Statistical and numerical contract

Let `Q` be the random-effect precision, `u_hat` the conditional mode, and
`eta = offset + X beta + Z u_hat`. The current evaluator maximizes, in `u`,

```text
conditional log likelihood(y | beta, u) - 0.5 u' Q u
```

and returns the first-order Laplace approximation

```text
logLik_L =
  log p(y | beta, u_hat)
  - 0.5 u_hat' Q u_hat
  + 0.5 log|Q|
  - 0.5 log|H_uu|

H_uu = Q + Z' W Z
```

The Gaussian normalization and Laplace integration constants cancel in this
expression. The Hessian is with respect to random effects, not the complete
fixed-plus-random coefficient vector. The likelihood is approximate marginal
ML, not exact integration or Gaussian REML.

For the currently assembled independent random terms, precision blocks are
`Sigma_term^-1 (Kronecker) A_term^-1`, using the relationship precision
supplied by model preparation. The outer optimizer estimates fixed effects
and the free descriptor coordinates.

The original foundation warm-started modes across objective evaluations;
milestone 1 uses deterministic zero starts for the outer objective. The native
evaluator still accepts an explicit mode start. Fixed-effect and covariance-
parameter standard errors are derived from a numerical Hessian over all
optimized parameters only when requested after fitting; they are withheld for
non-convergence, invalid curvature or detected covariance boundaries.

`re.form=NA` sets random effects to zero in predictions. For a nonlinear link,
this is **not** a population-marginal mean obtained by integrating over random
effects. Conditional mode curvature is likewise not automatically a complete
joint prediction-error covariance.

## 5. Validation record and limits of the evidence

Milestone 1 focused qualification:

- Expanded Laplace suite: **128 passing assertions**, without warnings or
  failures. Earlier combined Laplace/PQL/control/capability run: **370 passes**.
- Source-package build, including vignettes, succeeded.
- Final full-suite and CRAN-style check were stopped before completion when
  the user redirected work to M2. Incoming external DOI/URL checks had TLS
  connection resets. The full suite had reached late prediction/PXEM tests;
  two Laplace warnings were observed but no final warning details were
  obtained. This is not final-snapshot release qualification.
- PDF manual generation is blocked by unavailable `pdflatex`.

Milestone 2 final focused qualification:

- **282 passing Laplace assertions**; **394 passing assertions** combined with
  PQL, control/capability and Gaussian fixed-design regressions. No reported
  warnings, failures or skips on this final run. Legitimate unavailable-SE
  warnings in weak-curvature examples are explicitly expected and their
  diagnostic/output contracts are checked, not suppressed in fitting code.
- Independent R covariance-inverse/Kronecker objectives agree with native
  descriptor precision objectives to the tested `1e-10` tolerance; modes
  agree to `1e-9`. Finite-difference gradients cover both families, scalar
  and unstructured covariance, AR, crossed terms, random slopes and fixed
  covariance parameters.
- Regressions cover changed sparse patterns, identical-parameter/control
  mode reuse, fresh/cached agreement, strict `1e-12` mode-score convergence,
  binomial tails at predictors +/-40 and +/-700, serialization, public
  panel-budget parity and explicit insufficient-memory errors.
- Native compilation succeeded. The last C++ correction was compiled with
  `R CMD INSTALL --libs-only --no-test-load --no-configure --no-multiarch`,
  reusing existing Gaussian objects. Follow-up tests used `compile=FALSE`;
  no further full package recompilation was scheduled.
- The options-guide builder/checker, both modified Rd syntax checks and
  `git diff --check` passed. Gaussian engine bodies and the existing
  covariance/support headers have no source diff.
- Full package/CRAN-style qualification and cross-platform/catalog coverage
  are still not established by these focused results.

Post-build validation, from the package root:

```r
pkgload::load_all(".", compile=FALSE, quiet=TRUE)
testthat::test_dir(
  "tests/testthat",
  filter="mmes-laplace|mmes-pql|mme-control$|mme-capabilities|mmes-fixed-design",
  package="sommer", load_package="none", reporter="summary",
  stop_on_failure=TRUE
)
```

Recorded at the end of native implementation:

- Historical foundation Laplace suite: **42 passing assertions**, no warnings or
  failures.
- Focused existing PQL and capability tests passed.
- A dense relationship-precision regression converged after replacing the
  earlier L-BFGS-B outer optimizer with `nlminb`.
- The scalar analytic test checks the mode score, random Hessian determinant,
  and likelihood normalization/correction.
- Package options guide regenerated; whitespace checking passed.
- A full native-era suite run reported **8,612 passes, one failure, one skip,
  and no warnings**. Its failure was the new one-level factor prediction test
  against code loaded before the factor-level fix. The subsequent focused
  suite passed after the fix, including a missing local binding corrected in
  that rerun.

There is **not yet a clean full-suite run on the final source snapshot**, nor
a completed native-era `R CMD check`/cross-platform qualification. The earlier
clean full run for the discarded adapter is not native-engine validation.
Convergence of one relationship example, especially near a zero variance, is
not statistical calibration or scalability evidence.

Milestone 3 final focused qualification:

- **281 expanded-family assertions** and **282 retained M1/M2 Laplace
  assertions** passed, with no reported warnings on the final mathematical
  regression run. Including PQL (76), Gaussian fixed-design (24), multivariate
  PQL (21), control (5) and separately rechecked capability coverage (7),
  **696 focused assertions passed across the completed checks**.
- The combined run initially loaded capability metadata before its new
  empty-list configuration leaf was added. Its only remaining failure was
  that catalog mismatch; current metadata/leaf traversal was subsequently
  rechecked successfully in a focused seven-assertion run. Do not confuse
  this record with a clean full-package/CRAN-style check.
- All seven family/link cases match independent R density normalization,
  first scores and predictor curvature through its second derivative.
  Shape/dispersion mixed derivatives were checked independently. Full
  marginal scores use centered `1e-4` differences and an explicit maximum
  absolute discrepancy threshold of `2e-5`, including models with no free
  covariance parameters and mixed estimated nuisance parameters.
- Grouped-binomial `binm()` and two-column fits agree across all three
  links. Independent binomial normalization/Laplace correction, non-prefix
  filtering, trial retention, proportion-scale residuals and invalid-input
  contracts pass.
- Gaussian native likelihood agrees with independent exact marginal GLS
  (likelihood tolerance `1e-9`, coefficients `1e-7`), with known residual and
  random covariance.
- Seed 817, 300 observations/30 groups, generated NB2 shape 2.5 and Gamma
  dispersion 0.2 with random SD 0.25. Recovery regressions require shape
  discrepancy below 1.2, Gamma dispersion discrepancy below 0.07 and fixed
  effects within their declared tolerances. Estimated NB2 likelihood is
  at least the fixed-shape likelihood, to `1e-7`.
- Seed 207, 200 long-format NB2/Gaussian observations/20 groups, generated
  shape 2.5 and residual variance 0.7 with correlated latent covariance
  fixed at its known value. Fitted shape was about **2.09838** and Gaussian
  dispersion **0.696019**; nuisance/fixed full-score maximum absolute
  discrepancy at the optimum was **6.87394e-8**. These are deterministic
  recovery regressions, not population-wide estimator calibration.
- A single mixed fixture covers all seven links/families, positive trials,
  missing fixed covariates/trials, correlated random effects, fixed family
  parameters, row-family indexing, prediction and serialization. NB2 shape
  objects and Gaussian/Gamma residual variance/deviance conventions use
  fitted, not starting parameters.
- Native compilation succeeded using the existing Gaussian objects.
  The 92-option guide builder/checker, five modified/new Rd syntax checks,
  generated registration, whitespace checks and unchanged Gaussian/
  covariance-support/dependency-manifest checks passed.
- No complete package check, exhaustive covariance-catalog qualification,
  large-scale recovery study or cross-platform validation is claimed.

## 6. Gap assessment and ordered implementation milestones

### M1. Correctness and release hardening -- implemented, full qualification interrupted

The following checklist records the foundation gaps that motivated M1; the
execution update above describes their implemented resolution.

1. Propagate all supported preparation arguments consistently. The original
   bridge does not forward the supplied contrasts into `setupArguments`, and
   the dispatch does not carry the non-response missing-data policies into
   the Laplace preparation path. Store fitted contrast coding and training
   terms for reproducible new-data design construction.
2. Simplify and test observation/offset alignment, including missing
   fixed/random covariates as well as responses. The original concern about
   non-prefix omissions was based on misreading a logical mask as integer
   indices; it was not a demonstrated defect.
3. Correct native Poisson validation: the original C++ predicate combined
   negativity and non-integrality with AND rather than OR. The R entry point
   rejects invalid counts, but the native function must enforce its own
   contract too.
4. Return structured reasons for invalid mode/covariance evaluations. Audit
   broad covariance catches, the constant objective penalty, optimizer
   termination, separation, and mode iteration-limit behavior. Expected
   rejected proposals should be distinguishable from programming errors.
5. Diagnose Hessian definiteness and boundary fits explicitly. Do not silently
   turn negative variance estimates from an indefinite inverse Hessian into
   zero standard errors. Separate mode convergence, outer convergence,
   curvature validity, and boundary/singularity status.
6. Replace blanket descriptor-coordinate limits with descriptor-aware
   admissibility and boundary reporting. Verify warm-start/order invariance
   to the required likelihood and derivative tolerances.
7. Audit post-fit methods and controls: unsupported `predict()` aggregation,
   PEV, random-form specifications, and other options must not be silently
   ignored. Correct working-residual semantics and test all residual types.

**Acceptance:** targeted regression tests for each issue, independent scalar
and multivariate likelihood/Hessian checks, repeated-start agreement,
deterministic failure reporting, then a clean frozen-snapshot full suite and
`R CMD check`. Distinguish first-order approximation error from code error.

**Main-engine changes:** none required.

### M2. Modular native architecture and scaling seam -- core implemented and focused-validated

Implementation and focused validation are complete in the working tree.
Small profiling results are recorded below; large-scale qualification remains
separate from these results.

**Implemented seam**

- `.mmes_prepare_model()` uses `.mmes_fit(..., .prepareOnly=TRUE)` with the
  actual binomial/Poisson family, not the Gaussian/PQL inner-fit flag.
- A fitting-local external-pointer context owns cloned descriptors, designs,
  offsets, immutable relationship precisions/determinants and CHOLMOD state.
  It is never retained in the serialized fitted object.
- Descriptor covariance inverses/log determinants/derivatives use the
  existing `CovarianceFactorEngine` support without changing that support or
  either Gaussian engine. Identical descriptor coordinates reuse precision
  and derivative states. Complete block assembly is skipped if unchanged.
- CHOLMOD symbolic analysis is reused only when its sparse pattern matches.
  Determinant extraction copies the factor because `factor_to_sparse`
  consumes numeric factor storage needed for subsequent solves.
- Distinct proposals solve from zero. Only identical parameters **and** mode
  tolerance/iteration controls reuse a converged mode/score.
- Near likelihood rounding precision, a backtracking step must improve the
  mode score (or satisfy its requested tolerance), allowing only a
  machine-scaled likelihood discrepancy. This prevents strict derivative
  checks from failing on objective cancellation without relaxing mode
  tolerance or accepting an arbitrary mode.
- Exact trace calculations use dense solve panels of at most 32 columns,
  sized against `execution$denseMemoryMB`. This previously unused public
  Laplace control is now accepted and tested. This bounds panel
  workspace only: sparse Kronecker blocks, derivative blocks, family vectors
  and sparse factor fill still require memory. Large-scale matrix-free
  precision actions/selected inversion remain future optimizations.
- `nlminb` receives full marginal first scores; `optimHess` remains numerical
  for uncertainty. `laplaceDiagnostics$gradient` is the marginal log-likelihood
  score; `nativeWork` counts calls, mode evaluations, symbolic analyses,
  relationship factorizations and descriptor updates/reuses through the
  final optimizer score evaluation, before the uncertainty Hessian.

**First-score equations** (not Gaussian AI scores)

For `H = Q + Z' W Z`, `v = diag(Z H^-1 Z')`,
`a = 0.5 (dW/deta) * v`, and `c = H^-1 Z' a`:

```text
d logLik_L / d beta = X' [score - a + W * (Z c)]

d logLik_L / d theta =
  -0.5 u' D u + 0.5 d log|Q|/d theta
  -0.5 tr(H^-1 D) + c' D u
  where D = dQ/d theta
```

These include the dependence of the converged mode and curvature determinant
on every optimized parameter. Factor derivatives follow existing descriptor
support; a custom factor's own numerical derivative is not magically
converted into an analytic derivative.

**Small profiling evidence**

Four Poisson objective evaluations used 240 observations, seed 721 and a 1 MB
dense panel budget. There were 20 group levels with 12 records each, six
environment levels and 16 crossed levels cycled across rows. The covariate
cycled across 12 equally spaced values in [-1, 1]; counts were generated at
`exp(0.1 + 0.2*x)`. Fixed parameters were `(0.1, 0.2)` and covariance
parameters were the model's declared starting coordinates. Shape/slopes/AR
and sparse-relationship examples used `sigma2=0.5`.

| Model | Random coefficients | Parameters | Mode (s) | Extra native score (s) | Numerical score (s) | Max absolute score error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Crossed `group + crossed` | 36 | 4 | 0.010 | 0.002 | 0.034 | 1.15353e-8 |
| Random slopes `usm(leg(x, 1))` by group | 40 | 5 | 0.004 | 0.003 | 0.010 | 1.12192e-8 |
| AR environment covariance by group, rho=0.25 | 120 | 4 | 0.002 | 0.005 | 0.008 | 2.26564e-8 |
| Sparse relationship precision by group | 20 | 3 | 0.002 | 0.001 | 0.016 | 1.18853e-8 |

The sparse relationship was symmetric tridiagonal with diagonal 1.5 and
off-diagonal -0.2, marked as inverse input. Numerical scores used centered
`1e-4` steps and mode tolerance `1e-10`. All fresh/cached objective differences
were zero at recorded precision; fresh/cached scores passed `1e-9` agreement.
Each model used one symbolic analysis across the controlled derivative
proposals, and one invariant relationship factorization per term (two for
the crossed model). Asking for the score after the same objective reused the
mode rather than solving it again.

These single millisecond-scale measurements on a shared Linux/R 4.2.2
environment are descriptive, not portable speed-up or peak-memory claims.
For the 120-coefficient structured example, exact leverage/trace solves cost
more than the conditional mode: derivative solves/selected-inverse methods
are the demonstrated next scaling target. No changes to Gaussian AI are
needed to pursue it. These are not large genomic/pedigree benchmarks,
simulation recovery, or full-catalog/cross-platform qualification.

Original M2 scope/acceptance:

- Extract a preparation-only contract shared with Gaussian/PQL code instead
  of depending on the internal Gaussian working-model flag.
- Define reusable family output: log density, mean, score, curvature, and
  required higher derivatives. Separate family evaluation, mode solve,
  precision operations, and marginal objective.
- Use descriptor-native precision actions, log determinants, and derivatives
  where available, instead of repeatedly materializing covariance inverses
  and explicit Kronecker precision blocks in R.
- Cache invariant relationship factorizations/determinants; retain symbolic
  sparse analysis across objective calls when the pattern is unchanged.
- Add derivatives of the **Laplace objective**, including mode dependence and
  determinant correction. Validate against carefully controlled numerical
  derivatives before using them for optimization.
- Benchmark crossed random effects, random slopes, structured multi-level
  covariance, and sparse pedigree/genomic examples, with memory budgets.

**Acceptance:** objective/gradient agreement on independent test problems;
documented profiling showing the actual bottleneck and numerical parity
before promoting any optimization.

**Main-engine changes:** prefer extracting only demonstrated reusable
precision/factorization helpers into shared support. Any call-site edits to
`ai_mme_sp2()` must retain Gaussian objective/derivative parity. Do not route
GLMM updates through Gaussian AI scores. Keep `ai_reml_direct_sp2()` unchanged
unless a separately justified direct GLMM strategy is approved.

### M3. High-value response coverage -- agreed family set implemented and focused-validated

The user explicitly approved the full bounded family set below, rather than
only grouped binomial/NB2. Implementation is native and does not modify the
Gaussian coefficient engines or introduce a package dependency.

| Family | Link(s) | Response and nuisance semantics |
| --- | --- | --- |
| Binomial | logit, probit, complementary log-log | Binary or integer successes/trials; combinatorial constants retained. |
| Poisson | log | Non-negative integer counts; fixed dispersion one. |
| NB2 | log | Non-negative integer counts; `Var(Y)=mu+mu^2/theta`; positive shape estimated by default. |
| Gaussian | identity | Finite real response; positive conditional residual variance estimated by default; integration is exact. |
| Gamma | log | Strictly positive response; `Var(Y)=phi*mu^2`; positive dispersion estimated by default, shape `1/phi`. |

**Direct changes**

- `mme_glmm_family.h` uses small truncated Taylor kernels in the predictor
  (order four) and log nuisance parameter (order one). It supplies first
  scores, observed curvature, its first/second predictor derivatives and
  mixed nuisance derivatives without finite-differencing family scores.
  Stable logit tails, normal log-CDF/hazard recursion for probit and
  `expm1`/`log1p` complementary log-log calculations avoid naive probability
  subtraction. NB2/Gamma densities use R's scalar density/special-function
  routines for normalization, not another fitting package.
- Native objective parameters are ordered fixed effects, free random
  covariance coordinates, then free family log-shape/dispersion coordinates.
  Context caching includes every one of these parameters and mode controls.
  Family nuisance scores are evaluated even with **no** free random covariance
  coordinates.
- `nbinom2(theta=...)` provides the native R family. Existing negative-binomial
  family objects with the NB2 variance convention and log link are accepted,
  but are not required. NB2 shape and Gaussian/Gamma dispersion are estimated
  by default. `task$glmm$familyParameters` supplies explicit starting/fixed
  values; mixed traits use trait-named nested lists.
- `binm()` proportions/attributes and two-column binomial responses are
  converted to native success counts. `familym(..., trials="column")`
  explicitly supports success-count binomial rows among arbitrary other
  response types. Trials participate in row filtering and remain attached
  to retained rows. Binomial fitted values and residuals stay on the
  proportion scale.
- Long-format `familym()` assigns one native likelihood per observation.
  Correlated random covariance may connect traits; observation likelihoods
  remain conditionally independent. There is no PQL residual-reference-trait
  constraint. Structured residual correlation remains unsupported.
- Fitted family objects are updated to fitted NB2 shape. Summary/print,
  dispersion/nuisance tables, Pearson/working/deviance residuals, response
  prediction, serialization and AIC/BIC parameter counting are wired for
  all declared families. Approximate natural-scale nuisance SEs use the
  joint numerical Hessian and existing convergence/curvature/boundary gates.

For nuisance coordinate `lambda`, the full marginal first score is

```text
d logLik_L / d lambda =
  sum(d log p(y|eta)/d lambda)
  - 0.5 sum(v * dW/d lambda)
  - (Z c)' (d score/d lambda)
```

Here `v`, `c` and the mode are as defined in M2. The last term accounts for
the nuisance dependence of the conditional mode; it must not be omitted
when estimating NB2 shape or continuous-family dispersion.

**Validation:** expanded tests and directly related PQL/control/multivariate/
preparation regressions passed as recorded above. Small seeded recovery checks are
implementation regressions, not statistical calibration or assurance of
identifiability in every mixed-family/covariance model.

**Still outside this bounded M3 delivery:** beta-binomial, other links,
arbitrary custom native densities, dispersion regression formulas,
zero inflation/hurdles, structured residual correlation and broad
simulation/cross-platform qualification. M4 inference/prediction work and
M5 extended distributions remain separate.

Original recommended order:

Recommended order:

1. Grouped-binomial trials, with combinatorial constants and explicit
   trial/weight semantics.
2. Negative-binomial likelihood and shape estimation; then beta-binomial
   overdispersion if prioritized.
3. Additional links and continuous likelihoods where scientifically useful,
   with family-specific derivative and tail-stability tests.
4. Mixed-family `familym()` integration using observation-family indexing and
   a coherent joint latent-effect model.

**Acceptance:** independent density/derivative checks, known-parameter
simulation recovery, offsets/trials/missingness tests, and declared
identifiability limits. Existing PQL negative-binomial support is not native
Laplace negative-binomial support.

**Main-engine changes:** none expected in either Gaussian objective.

### M4. Inference, prediction, and simulation -- implementation in place; broader qualification outstanding

- Expose validated fixed-effect covariance and natural-scale
  covariance-parameter uncertainty, with boundary/curvature warnings.
- Add conditional new-data prediction, explicit unseen-level handling, and
  separate conditional, zero-random-effect, and population-marginal means.
- Define joint fixed/random uncertainty and coefficient PEV approximation;
  do not present `H_uu^-1` alone as full prediction-error uncertainty.
- Add Wald contrasts and carefully qualified nested-model comparisons.
  Variance-component boundary LRTs require nonstandard calibration or
  parametric bootstrap, not an unconditional chi-square reference.
- Implement native simulation and randomized quantile residuals as tools for
  model checks; dedicated diagnostics for dispersion, excess zeros, residual
  dependence, and broad predictive calibration remain to be qualified.

**Acceptance:** simulation-based coverage/calibration and post-fit method
compatibility tests, not merely finite standard errors.

**Main-engine changes:** no Gaussian objective changes; reuse prediction
helpers only after documenting compatible covariance semantics.

**Execution update:** the joint-information, deferred uncertainty, prediction
and sparse simulation tranches are implemented in the working tree. This is
not full release qualification.

- `vcov()`, `summary()`, and `covparams_mmes_se()` calculate numerical-
  Hessian uncertainty after fitting. `mmes()` no longer computes or stores
  the optimizer Hessian during native Laplace fitting. PEV requests during
  fitting are rejected unless `output$inverse=0`; `postPEV()` performs the
  sparse joint fixed/random information solve after fitting. The LMM output
  controls remain unchanged.
- Free covariance-coordinate uncertainty is transformed to natural-scale
  standard errors in post-fit summaries; fixed parameters remain explicitly
  marked. All parameter uncertainty stays behind convergence, positive-
  curvature and boundary gates.
- `anova(fit)` provides term-wise asymptotic Wald chi-square tests.
  `anova(smaller, larger)` is limited to strictly nested fixed effects with
  identical response, offset, family/trial, random design, relationship
  precision and covariance structure. Boundary fits and variance-structure
  comparisons are refused; no ordinary chi-square variance-component LRT is
  exposed.
- New-data predictions distinguish conditional fitted modes
  (`re.form=NULL`) and zero-random-effect means (`re.form=NA`). Conditional
  prediction reports unseen levels and requires known levels. Population
  response means (`re.form="population"`) use 24-point Gauss-Hermite
  integration only for independent random-intercept terms with identity
  relationship precision; unsupported covariance structures error explicitly.
- `predict(..., se.fit=TRUE)` reports post-fit conditional/zero-random-effect
  uncertainty and restricted population-mean uncertainty. `cov.fit=TRUE`
  additionally returns full prediction covariance under the configured
  memory budget. The approximation conditions on fitted covariance and
  family parameters for conditional predictions; population predictions
  propagate the joint fixed/covariance-parameter curvature using numerical
  prediction derivatives. Both routes use `execution$denseMemoryMB` and
  `output$panelSize` where applicable.
- `simulate()` draws family-specific response replicates conditional on the
  fitted random modes by default, or draws new Gaussian random effects from
  fitted covariance blocks with `re.form="population"` using sparse
  relationship-precision factors. Estimated workspace and returned
  simulations are guarded by `maxMemoryMB`, defaulting to
  `execution$denseMemoryMB`. Draws include grouped-binomial proportions and
  fitted NB2/Gaussian/Gamma nuisance values. `residuals(type="quantile")` adds
  randomized probability-integral-transform residuals for discrete responses
  and fitted-CDF residuals for Gaussian/Gamma. These are conditional response
  checks or model-based random-effect simulations, not a parametric bootstrap
  or population predictive calibration.
- `postPEV()` computes approximate Laplace fixed/random PEV output only after
  fitting: mode 1 computes random-effect
  prediction-error diagonals using bounded sparse solve panels, and mode 2
  additionally returns the full joint fixed/random conditional
  observed-information inverse. The random block is the conditional PEV
  approximation with fitted covariance and family parameters held fixed;
  mode 2 is guarded by `execution$denseMemoryMB`. Post-fit calculations
  do not mutate likelihood inputs, invoke either Gaussian coefficient engine,
  or require a package dependency.
- Conditional new-data preparation exposed a one-row sparse residual design
  edge case; `vsm()` now normalizes the one-column design to sparse general
  form before reading its triplets. This does not change either Gaussian
  coefficient objective.
- Focused native Laplace and expanded-family regression files pass on the
  current inference/simulation/prediction changes. Tests cover guarded
  `vcov()`, term-wise `anova()`/`wald_mmes()`, strictly nested fixed-effect
  likelihood comparisons before and after deferred PEV calculation, deferred
  PEV modes, conditional/zero/population prediction covariance, one-row
  prediction, named unseen-level errors, seed-preserving sparse prior
  simulation, grouped mixed-family draws, and randomized quantile residuals.
  An exact Gaussian marginal-GLS comparison verifies `vcov()`; a seeded
  1,000-replicate known-covariance simulation checks nominal fixed-effect
  interval coverage. A seeded population-response simulation also agrees
  with the restricted integrated mean. The restricted
  population mean agrees with independent numerical integration for a
  seeded random-intercept case. Rd parsing and whitespace checks pass.
- Final focused Laplace/PQL/control/capability/fixed-design run passed **765
  assertions**, with one optimizer warning for false convergence in a
  structured-covariance regression. The full testthat suite reached
  completion but failed in the unchanged `test-mme-threads.R` cleanup:
  `Sys.setenv()` received unnamed arguments. This failure is outside the GLMM
  changes; the full suite is therefore not clean.
- `R CMD build --no-build-vignettes .`, the options-guide checker (92 options,
  17 control boxes) and individual Rd parsing passed. Source-directory
  `R CMD check --no-manual --no-build-vignettes .` stops at package metadata
  because `DESCRIPTION` lacks the required explicit `Author` field; a
  completed package check is not claimed.
- General population integration beyond independent random intercepts,
  variance-component parametric-bootstrap calibration, broad non-Gaussian
  interval-coverage studies, and dedicated excess-zero/residual-dependence
  diagnostics remain outstanding. Full package checks and cross-platform
  qualification are also outstanding; this bounded M4 implementation is not
  a claim of full release qualification.

### M5. Broader competitive feature roadmap

- **Initial zero-inflation slice implemented:** `zipoisson()` supplies a
  Poisson log-mean predictor plus an independently estimated/fixed
  intercept-only logit structural-zero probability. The native likelihood,
  parameter score/curvature, `familym()`, response-scale means, simulation and
  randomized quantile residual CDF are connected. The zero probability is
  included in the family-parameter summary, post-fit uncertainty and AIC/BIC
  parameter count. This is native-Laplace-only and does not add a package
  dependency. It is not a hurdle model, zero-inflated NB2, dispersion model,
  or general two-formula model. The zero-component intercept can be
  confounded with the count intercept and requires further identifiability
  and calibration study.
- Focused native Laplace, expanded-family, and capability regressions pass
  with no failures. The expanded-family suite checks independent zero-
  inflated densities, predictor/parameter derivatives and full marginal
  gradients; end-to-end single- and mixed-family fits; natural-scale
  uncertainty; fitted and predicted means; seeded simulation zero-frequency
  agreement; quantile residuals; PQL rejection; and fixed `familyParameters`.
  One existing structured-covariance regression warned that its optimizer
  reported false convergence; it is unrelated to the ZIP tests. These are
  implementation regressions, not broad statistical calibration.
- Initial reproducible simulation calibration completed **560 fits** across
  16 settings: 20 versus 60 groups (six observations/group), underlying
  Poisson rates 0.5 versus 2, fixed random-intercept variances 0.0001 versus
  0.16, and true zero inflation 0 versus 0.25. Each zero-inflation-null
  setting had 20 replicates; each positive-inflation setting had 50. All
  optimizer fits converged. This conditions on the supplied, known random
  variance; it does not assess variance-component estimation.
  - Under true zero inflation of 0.25, mean estimates ranged from 0.211 to
    0.246. Approximate 95% delta/Wald coverage for the zero probability was
    0.90 to 1.00 among fits with available curvature; fixed-effect slope
    coverage was 0.88 to 0.98. Some sparse/low-rate settings had 16% of
    zero-probability SEs unavailable because of boundary or curvature gates.
  - Under true zero inflation of 0, mean estimates ranged from 0.0039 to
    0.043, with RMSE 0.010 to 0.089. The zero-probability boundary was flagged
    in 45% to 75% of fits; valid post-fit zero-probability SEs were available
    in only 20% to 55%. These are expected nonregular boundary cases, but also
    show that optimizer convergence alone is not evidence that the extra
    component is identifiable.
  - The experiment is a screening study, not final coverage evidence: only
    20/50 replicates per setting, and the zero-inflation-null boundary and
    selected finite-curvature subsets make ordinary Wald coverage incomplete.
- The zero-inflation-null limit is also regression-tested against the
  Poisson likelihood, coefficients and fitted means with a fixed probability
  of 1e-8. Existing mixed-family regressions exercise ZIP alongside the other
  supported families, including fixed per-trait `familyParameters`.
- **Expanded baseline evidence:** an external, frozen-source study completed
  3,900 datasets, 11,700 primary fits and 780 alternative-start ZIP fits.
  Correct constant-pi ZIP at moderate count rate recovered mean pi 0.252
  (truth 0.25), with 94.7% interval coverage and 100% uncertainty availability.
  Sparse low-rate uncertainty availability was only 55.3%; exact Poisson
  cases remained nonregular. Covariate-dependent inflation produced count
  slope biases -0.058 and -0.149 under the constant-pi model. This evidence
  justified a bounded fixed-effects zero-formula prototype, not zero random
  effects or universal mixture identifiability.
- **Bounded formula slice implemented:** `zipoisson(zero=~covariates)` uses
  a separate fixed-effects logit design and optional offset, jointly fitted
  with the count and covariance parameters. Default `zero=~1` retains the
  scalar representation. Trait-aware missingness, factor/transform-aware
  prediction, likelihood comparison signatures, simulation, residuals,
  parameter counting and joint post-fit uncertainty are connected.
  Coefficients appear in `zeroCoefficients` on the logit scale.
  PEV and all uncertainty remain post-fit; no package dependencies were added.
- Formula qualification uses an externally registered contract and 1,800
  paired datasets (300 per scenario), with coefficient recovery, coverage,
  availability, constant-pi false-positive behavior, population function
  recovery and 360 alternative-start fits, for 3,960 fits in total.
  All source hashes and 1,800 legacy log-likelihood reproductions passed
  verification (maximum likelihood difference zero). Independent checks
  recomputed coverage denominators, Wilson intervals, paired Monte Carlo
  moments and population references.
  - **36/41 registered gates passed.** The three moderate-rate scenarios
    (rate 2, 60 groups, six records/group, zero slopes 0/0.3/0.6) passed every
    gate: 900/900 converged fits and joint SEs. The constant-pi zero-slope
    false-positive rate was 15/300 = 5.0%. At true zero slope 0.6,
    count-slope bias improved from -0.058 to +0.004.
  - **Low-rate qualification failed:** with rate 0.5 and zero slope 0.6,
    convergence was 293/300 = 97.7% (required 98%), joint SE availability
    was 244/300 = 81.3% (required 95%), zero-intercept bias was -0.412 and
    zero-slope bias +0.153 (required absolute bias <=0.10). These are five
    failed gates, not a reason to loosen the contract.
  - Across primary and sensitivity fits, 67 explicit numerical-boundary
    errors occurred in stress/low-rate scenarios; none occurred in the three
    moderate-rate scenarios. Failures and unavailable intervals remain in
    their planned denominators. Ordinary Wald inference at exact Poisson
    and sparse-mixture boundaries remains unqualified.
  - Focused native-family/capability regressions and formula/shared Gaussian
    preparation/prediction regressions pass; the existing structured-
    covariance false-convergence warning remains. Rd checks, the options
    catalog, selected article renders and source build pass. No full package
    check or cross-platform release qualification is claimed.
  - **Decision:** retain the explicit, bounded fixed-effects formula slice
    with a moderate-rate calibration envelope and clear low-count warnings.
    Prioritize identification/boundary diagnostics and calibrated inference
    before adding mixture components. The complete external evidence lives
    in the sibling `sommer-glmm-m5-validation/formula-v2` study directory,
    not in the package source tree.
- Keep zero random effects, hurdle components and zero-inflated NB2 deferred.
  A successful formula optimizer is not evidence that ordinary Wald
  uncertainty is valid in sparse, separated or boundary mixtures.
- Model distribution-specific dispersion/shape using a separate formula;
  this is not a Gaussian residual variance model relabeled for counts.
- Add component-specific random effects only after the basic multi-component
  likelihood and identifiability diagnostics are validated.
- Consider adaptive Gaussian quadrature for eligible low-dimensional or
  independent-group integrals as an accuracy check/alternative, not a generic
  replacement for sparse high-dimensional Laplace.
- Treat correlated non-Gaussian observations/latent residual processes as a
  separate model-design project. An arbitrary Gaussian `rcov` or precision
  `W` does not automatically define a valid count/binomial joint likelihood.
- Revisit direct formulation, accelerations, and very-large-model approximate
  traces only when user demand and measured bottlenecks justify them.

These are broader competitive capabilities, not prerequisites for an honest,
well-tested initial native binomial/Poisson release.

## 7. Recommended next deliverable

### Staged ASReml-gap roadmap

The approved sequence is reliability, ordered categorical responses, then
nominal multinomial responses. These are distinct milestones, not aliases
for one scalar family extension.

1. **Reliability foundation implemented:** `boot_mmes()` provides opt-in
   post-fit population parametric bootstrap for parameter-recovery summaries
   and nested count-fixed-effect likelihood-ratio tests. It reuses stored
   native likelihood designs/specifications, including mixed families,
   grouped binomial trials and ZIP zero formulas. Refits use native first
   derivatives without numerical Hessians or PEV. Covariance/family draws
   retain their declared optimizer coordinates rather than silently
   masquerading as natural parameters.
   - Simulations draw new random effects. Sparse relationship factors are
     reused across replicates, and their caches no longer modify fitted
     relationship metadata.
   - Every planned replicate has diagnostics for convergence, boundaries,
     likelihoods, scores, rejected evaluations, warnings and errors.
     Failed refits are not replaced. Any failure makes the ordinary p-value
     unavailable; failure-sensitivity bounds retain the planned denominator.
   - No Gaussian fitting controls, LMM routes, dependencies or fit-time
     uncertainty behavior changed. The bootstrap's outer limit is a
     post-fit argument; inner mode/resource settings come from the fit.
   - This does not yet provide covariance-component or structural-zero
     selection tests, profile intervals, automatic separation correction,
     or release-level calibration. A small external null-behavior screen
     is descriptive evidence, not proof of 95% coverage.
   - The external `bootstrap-v1` screen completed 40 Poisson and 40 grouped-
     binomial null datasets, with 49 paired replicates each: 3,920 pairs
     (7,840 native refits), all successful. Rejection at 5% was 0/40 for
     Poisson and 1/40 for grouped binomial; Wilson intervals were
     [0, 0.088] and [0.004, 0.129]. The 0.02 p-value resolution and small
     dataset count preclude a general calibration claim. Algorithm hashes
     and every plus-one p-value/denominator passed verification.
   - Bootstrap/ZIP, native-family and shared-preparation regressions pass;
     the pre-existing structured-covariance false-convergence warning
     remains. Rd checks, selected article renders and source build pass.
2. **Ordered categorical public implementation available experimentally:** cumulative logit/probit with
   a single predictor, explicitly ordered responses, identified location,
   ordered thresholds, category-aware prediction/simulation and joint
   post-fit uncertainty. Implementation gates pass; broader recovery/coverage
   calibration remains a qualification gate, not a claim implied by fitting.
   - **First likelihood slice implemented, internal only:** the native family
     evaluator accepts fixed, finite, strictly increasing cutpoints and
     integer categories `1,...,K`, with
     `P(Y <= k | u) = F(cutpoint[k] - eta)`. Positive predictor coefficients
     favor higher categories. The link scale is fixed to standard logistic
     or standard normal; no ordinal residual dispersion is estimated.
     It returns log category densities, expected category codes, predictor
     scores, observed curvature and its first two predictor derivatives.
     Logistic intervals use a stable factored identity; probit intervals use
     log-tail differences and a bounded narrow-interval density expansion.
     No numerical Hessian, PEV, dependency or Gaussian-path change is needed.
   - **Joint cutpoint estimation implemented internally:** free ordinal
     specifications now contribute one coordinate per cutpoint rather than
     using the one-scalar dispersion slot. The first coordinate is the first
     cutpoint; subsequent coordinates are log-positive successive gaps.
     Native updates reject overflow, underflow and gaps lost to floating-point
     rounding. Analytic endpoint density/probability ratios supply cutpoint
     scores and mixed predictor/curvature derivatives. Marginal scores include
     random-mode and log-determinant dependence, and transform to log-gap
     coordinates by the chain rule. A common cutpoint shift uses the negative
     predictor-shift identity to avoid narrow-interval cancellation.
     Per-family cutpoint vectors are aggregated without an observation-by-
     cutpoint derivative matrix. At this development stage public
     `mmes(family=ordinal())` was kept disabled. The constructor's
     variance/deviance/GLM methods now error explicitly instead of pretending
     ordinal data have Bernoulli variance or zero deviance.
   - Native objective construction checks location identification once, not
     inside optimizer evaluations. It rejects any linear combination of fixed
     columns confounded with free ordinal trait-location indicators, including
     equivalent intercepts and cross-trait indicator combinations. The dense
     rank check is admitted through the existing memory budget. Every declared
     category must be observed when estimating cutpoints; fixed-cutpoint
     evaluation does not impose that restriction. Other native families retain
     their existing parameter ordering and score machinery.
   - **Response/control preparation implemented internally:** the common
     preparation entry point recognizes a single ordinal family and reuses
     the existing observation map, missing-data policies, random structures
     and fixed design. Ordered labels/declared levels are preserved; plain
     numeric codes must be contiguous integers `1,...,K`. Unordered factors,
     characters, matrices and invalid codes error before row filtering.
     Every category must remain observed for estimated cutpoints; fixed
     cutpoints can retain absent ordered-factor levels.
   - Internal preparation uses existing `task$glmm$familyParameters` keys
     `thresholds` and `estimate`, with control thresholds taking precedence
     over constructor thresholds. Estimated cutpoints default to empirical
     cumulative-link quantiles computed after filtering. Fixed cutpoints
     require explicit values. Exactly `K-1` cutpoints are required. Estimated
     cutpoints remove only the existing intercept column (not formula `-1`,
     which would change factor contrasts). Remaining constant confounding and
     fixed-design rank deficiency error. Fixed cutpoints retain an intercept.
     Column selections, partitions, original formula/labels, contrasts and
     levels are retained for subsequent public prediction integration.
     Ordinal proposal updates and parameter counts use the same gap
     coordinates as the native objective.
   - Preparation initially kept public fit assembly gated until all downstream
     interfaces could be connected. Ordinal `familym()` remains explicitly
     unsupported.
   - Preparation-to-optimizer testing exposed an exactly-zero-curvature sparse
     iterator invalidation in native GLMM mode-Hessian scaling. Replaced both
     iterative and final native row-scaling loops with sparse diagonal
     multiplication, without changing Gaussian engines or mode objectives.
     Saturated ordinal and binomial regressions require finite valid results,
     not exceptions or suppressed errors.
   - **Public surfaces connected:** single-family cumulative-logit/probit
     fitting now assembles complete objects. Fitted/predicted distributions
     default to labelled category probability matrices; latent predictor and
     coding-dependent expected code are explicit types. Simulation preserves
     ordered labels, residuals use the category distribution or explicitly
     defined numeric-code moments, and likelihood counts include every free
     cutpoint. Comparison signatures preserve category labels and fixed
     thresholds; bootstrap refits retain gap coordinates and report missing-
     category failures. Mixed-family ordinal fitting remains unsupported.
   - Post-fit joint marginal-Hessian uncertainty transforms gap coordinates
     to natural cutpoint SEs. Prediction propagates all parameter cross-
     covariances and mode dependence, adding conditional random-mode variance
     when applicable. Probability SEs retain matrix shape; full covariance
     is column-major with category/row labels. Ordinal postPEV uses the same
     joint-parameter/mode decomposition instead of conditioning on estimated
     cutpoints. No fit-time numerical Hessian/PEV or R dependency is added.
   - Population probabilities use exact probit normal convolution and the
     existing 24-node Gaussian quadrature for logit, within existing random-
     covariance eligibility. Independent integration, probability, binary-
     limit and `MASS::polr()` likelihood tests cover the public interfaces.
   - **Acceptance gates:** independent category probabilities must sum to one
     within `1e-12`; binary-category limits must match native binomial
     likelihoods/derivatives within `1e-12`; predictor/mixed derivatives must
     match numerical references within `1e-6`, and full marginal gradients
     within `2e-5` in interior fixtures. Include both tails, narrow cutpoints,
     invalid responses and cutpoints. Joint-estimation integration must also
     establish identification, independent fixed-effect likelihood agreement,
     recovery, sparse-category failure reporting and external coverage gates
     before claiming a supported ordinal GLMM. Kernel checks alone do not
     qualify inference or fitting.
   - **Initial validation:** an isolated native package build and focused
     ordinal-kernel, native-Laplace, expanded-family, bootstrap and Gaussian
     fixed-design suites completed with 1,023 passing assertions, no failures
     or skips, and the known structured-covariance false-convergence warning
     in `test-mmes-laplace.R`. Both modified Rd files parse and check cleanly.
     Tail/narrow-probit derivatives also agree with independent integrated
     normal-density derivative references. These are implementation checks,
     not statistical calibration or a full cross-platform package check.
   - **Joint-estimation validation:** an isolated build and the same focused
     suites now complete with 1,195 passing assertions, no failures or skips,
     and only the known structured-covariance warning. Full marginal-gradient
     differences are checked against the absolute `2e-5` gate for both links,
     binary/multiple cutpoints, gaps down to `1e-5`, and multiple ordinal
     families alongside estimated Gaussian dispersion. Optimized likelihoods
     agree with independent `MASS::polr()` fits in a fixed vanishing-random-
     variance limit within `1e-5`, with coefficients/cutpoints within `1e-4`.
     Tests also cover equivalent intercepts, joint trait-location confounding,
     no fixed columns, memory admission, missing categories, parameter-vector
     lengths and recovery after rejected proposals. These do not establish
     random-variance recovery, ordinal interval coverage or public fit support.
   - **Preparation validation:** isolated installation and the expanded focused
     suites complete with 1,318 passing assertions, no failures or skips,
     and the same known structured-covariance warning. Ordered labels,
     retained-row initialization, fixed/estimated controls, missing levels,
     contrasts, formula environments, offsets, intercept-only designs and
     invalid schemas are checked. Prepared factor/offset models optimized
     through the native objective agree with independent `MASS::polr()` fits
     in the vanishing-random-variance limit at the existing `1e-5` likelihood
     and `1e-4` coefficient/cutpoint gates. The new Rd content parses and checks.
     This preparation-stage validation did not claim public fit or
     statistical-coverage qualification.
   - **Public-interface validation:** isolated installation and ordinal kernel,
     preparation/public API, native Laplace/families, ZIP formula, bootstrap,
     PQL, fixed-design and Gaussian ML/REML regressions complete with 1,720
     passing assertions, zero failures/skips, and only the known structured-
     covariance false-convergence warning. Public ordinal tests include both
     links, independent `polr()` likelihood limits, binary equivalence,
     fixed/estimated/intercept-only fits, filtered category schemas, new-data
     factor/offset coding, known/unseen random levels, simulation frequencies
     and RNG restoration, population integration, independent joint delta-
     method covariance, cutpoint SEs, PEV, bootstrap and comparison eligibility.
     Mocked uncertainty routines confirm fitting never invokes post-fit
     parameter uncertainty, PEV or `optimHess()`. Rd checks, selected article
     renders and generated option-guide checks pass. No Gaussian native
     solver policy or dependency was changed. Broader ordinal recovery/
     coverage calibration remains pending; these are implementation tests.
   - **Targeted ordinal calibration completed (external `ordinal-v1` study):**
     frozen source, 3,000 paired datasets (300 per link/design combination),
     6,900 fits including higher-order `ordinal::clmm(nAGQ=9)`, AGQ1 and
     alternate-start sensitivities. The external reference was already
     installed; no package dependency was added. Four fresh-process workers
     avoided fork/BLAS locks. The contract, source/binary hashes, deterministic
     checkpoints, raw fit records, CSV summaries, independent verifier and
     HTML report live outside the sommer source tree.
   - Both links were tested with four categories and random intercepts:
     60 groups x 8 records and 100 x 10 regular designs, plus short-cluster,
     sparse-category and near-zero-variance stress designs. All 1,200 regular
     primary fits converged; joint SE availability was 99.7--100%. All
     registered fixed-effect/cutpoint bias and natural random-variance bias
     gates passed. Regular alternate starts produced no likelihood gain
     exceeding `1e-4`. All numerical Hessian/interval work remained post-fit.
   - **The full targeted envelope did not qualify:** 18 of 712 registered
     regular gates failed, without changing thresholds after observing data.
     Fourteen concern probability calibration (three RMSE cells above 0.04
     and eleven Wilson-coverage checks), three indicate conservative parameter
     coverage, and one concerns AGQ1 reference agreement. Regular pointwise
     probability coverage ranged from 90.3% to 98.0%; every raw coverage-range
     gate passed, but some Wilson intervals excluded nominal 95%.
     These many correlated pointwise gates are not a simultaneous-coverage
     test; failures still block an unqualified validation claim.
   - A separate post-study scalar random-intercept Laplace calculation at
     each implementation's own estimates checked 120 regular sensitivity
     datasets. Sommer's maximum absolute likelihood difference was below
     `7.4e-9`; the AGQ1 reference differed by up to `1.76e-4`, including one
     reference nonconvergence. This localizes the failed reference-comparison
     gate to reference precision rather than demonstrating a native likelihood
     bug. The original gate remains failed; no record or tolerance was replaced.
   - Sparse fit availability was 78% (logit) and 14% (probit), exactly matching
     66/258 datasets with an absent declared category and intentional
     identification rejection. Every primary dataset with all categories
     observed converged. Joint SE availability across all planned sparse
     datasets was only 51.7%/9.7%; weak-variance designs supplied SEs in only
     52.3%/54.3% despite 100% optimizer convergence. Short-cluster logit
     variance recovery also differed from higher-order quadrature. These
     regimes are not qualified for ordinary interior inference.
   - Independent audit verified every dataset/fit count, source/binary hash,
     seed/checkpoint reproducibility, probability truth and sum-to-one,
     parameter/probability coverage numerator and denominator, and registered
     gate decision. The study does not qualify random-effect PEV coverage,
     fixed-cutpoint inference, general covariance, or cross-platform release.
     The next ordinal statistical objective is bounded post-fit probability-
     interval calibration (current Wald versus transformed/bootstrap intervals),
     with weak-information diagnostics, before broad inferential claims or
     a nominal-multinomial expansion. Do not silently change the Laplace
     estimator or drop failed/missing-category replicates.
   - **Held-out probability-interval calibration completed:** an independent
     frozen-source run generated 1,200 fresh datasets (300 per regular scenario)
     with seed base `328000000`, distinct from the original study. The same
     fitted probabilities/SEs were used for paired 95% Wald, logit, probit and
     arcsine delta intervals. All fits converged; SE availability was again
     99.7--100%. Across 96 registered category/grid cells per method, failures
     were 8 (Wald), 3 (logit), 4 (probit), and 5 (arcsine). No method passed
     the full envelope. Rich scenarios passed all cell gates, but the
     lower-information regular designs did not. The logit transformation
     reduced gate failures; that is not a universal coverage correction or
     grounds for selecting a new default after seeing results.
   - The independent interval verifier checks all 115,200 stored intervals
     against raw fits, transform endpoint formulae, coverage numerators,
     denominators, Wilson intervals, registered decisions and hashes.
     Transformed intervals are pointwise local approximations, not bootstrap
     intervals, joint-simplex regions or future-response intervals. No public
     uncertainty API/SE calculation changed. A genuine refit-bootstrap interval
     study remains pending. External contract, runner, audit and HTML report
     reside in the sibling validation project under `ordinal-interval-v1`.
3. **Nominal multinomial milestone started internally:** joint category-specific
   predictors, reference-category constraints and observation-block
   curvature. Independent binomial `familym()` rows are not a substitute
   for this joint likelihood. Extend the native block derivative contract
   before exposing a nominal family.
   - The first concrete slice is an isolated categorical softmax kernel:
     `K-1` predictors with the last category as zero-reference, one joint
     observation likelihood, full category-block score/curvature and the
     first two directional derivatives of that curvature. These are third/
     fourth categorical cumulant contractions, avoiding category-fourth-order
     tensor storage. Each output block uses predictor-category coordinates,
     not independent Bernoulli rows. Category permutation requires rebasing
     predictors to the new reference; adding an unidentifiable reference
     intercept is not permitted by this representation.
   - The internal evaluator is separate from scalar GLMM and Gaussian engines.
     It accepts finite `n x (K-1)` predictor/direction matrices and categorical
     codes `1,...,K`, uses stable log-sum-exp, and validates dimensions,
     response codes, numerical range and dense workspace before allocations.
     It returns `n x K` probabilities, `n x (K-1)` scores, and three
     `(K-1) x (K-1) x n` curvature arrays; the directional sign convention
     is for *negative* log-likelihood curvature. No fitting, Hessian inversion,
     PEV, dispersion or public nominal-family constructor is exposed.
   - **Conditional-mode slice implemented internally:** a tagged persistent
     context accepts `n x (K-1)` fixed predictors, category-major stacked
     sparse `n(K-1) x q` random design and supplied fixed SPD precision `Q`.
     It computes the joint score `Z' vec(score) - Q u`, the full observation-
     block Hessian `Q + Z' W Z`, Newton ascent with step halving, and the
     normalized-prior Laplace likelihood. This is not a marginal parameter
     optimizer. A structural Hessian envelope retains numerical zeros across
     saturation/cancellation, preserving CHOLMOD symbolic analysis across
     iterations and repeated predictor proposals. Workspace/index admission,
     tagged-pointer guards, explicit numerical failure reasons and work
     counters are included. The workspace budget covers block and structural
     assembly estimates, not an absolute bound on sparse-factor fill-in.
     No Hessian inverse, PEV or other uncertainty is computed during this work.
   - **Analytic marginal parameter scores implemented internally:** the new
     evaluator accepts a proposed SPD precision, sparse category-major
     fixed-predictor Jacobian and a list of first precision derivatives on
     caller-supplied optimizer coordinates. It re-solves the mode and returns
     log-likelihood scores in fixed-Jacobian-column then precision-derivative
     order. The score includes both normalized-prior determinant derivatives
     and mode-dependent observation curvature in the Hessian determinant.
     Small observation-wise third-cumulant contractions and one adjoint solve
     avoid per-parameter mode Jacobians. Precision traces and observation
     contractions use bounded sparse-solve panels, not full inverses, PEV or
     fitted-parameter uncertainty. Fourth cumulants are not needed for first
     marginal scores and remain available for later second-derivative work.
     Both factors reuse symbolic analysis for unchanged structure and refresh
     it for actual precision sparsity changes. Disabled gradients perform no
     score solves; failed modes return no gradient. Derivative dimension,
     symmetry, finiteness and workspace guards raise explicit errors.
   - **Nominal response/design preparation implemented internally:** unordered
     factor labels or contiguous plain numeric codes are explicitly mapped to
     `K-1` category contrasts and a last-position zero-reference. The reference
     defaults to the last declared level or is selected by label; original
     level order is retained in metadata. Ordered/character/matrix responses,
     trials and missing retained categories fail explicitly. The initial fixed
     design gives each non-reference category its own copy of the complete
     formula basis, retaining intercepts, transformations, contrasts and
     training levels. Formula offsets are repeated across contrasts. Other
     fixed-sharing and category-specific offset formulas remain outside this
     slice. A category-major observation map and common missing-data filtering
     prevent treating category rows as independent observed responses.
   - Random designs reuse existing formula/descriptor preparation on an
     internal stacked data frame exposing `.nominalCategory`. Plain `~group`
     explicitly means a shared adjustment across non-reference contrasts;
     `vsm(ism(.nominalCategory), ism(group))` gives independent category-specific
     adjustments with a common scale, and `usm(.nominalCategory)` permits
     category covariance. No random formula is silently rewritten. Existing
     relationship precision alignment and fixed covariance flags are retained.
     This bounded internal preparation currently requires an explicit data
     frame. It returns preparation metadata, not an `mmes` fitted object.
   - **Covariance-coordinate plumbing implemented internally:** a tagged
     nominal objective reuses `LaplacePrecision` and its descriptor derivatives
     for fixed coefficients followed by free covariance coordinates, excluding
     the residual placeholder. First derivatives retain their term-block
     offsets and analytic prior log-determinant derivatives rather than
     expanding every term to a full `q x q` derivative matrix. Repeated
     likelihood/score calls at identical parameters and mode controls reuse a
     successful mode/score; failed evaluations are never cached as successes.
     Caller-owned descriptors/results are not mutated. The existing modular
     `mmes_options`/`mmes_control()` validation supplies Henderson Laplace ML,
     mode budgets/tolerance, dense workspace and BLAS-thread controls.
     Inverse requests, unsupported family settings and unused Gaussian
     controls fail explicitly. Scalar GLMM/Gaussian routes and policies remain
     untouched; no dependencies or fitting-time uncertainty are added.
   - Next implementation gates: complete nominal fit/optimizer assembly and
     all downstream outputs, prediction, simulation, model comparisons,
     bootstrap and strictly post-fit uncertainty before exposing a public
     nominal family. Statistical calibration remains a separate release gate.
     Ordinal interval findings neither qualify nor invalidate this distinct
     joint nominal likelihood.
   - **Initial multinomial kernel validation:** isolated source installation and
     focused multinomial/ordinal kernel, preparation/public API, native Laplace/
     families, ZIP, bootstrap, PQL and Gaussian fixed-design/ML/REML regressions
     pass 1,887 assertions with zero failures/skips and the same known structured-
     covariance false-convergence warning. Tests cover independent probabilities,
     joint cross-category curvature, binary equivalence through order four,
     directional and mixed fourth-order finite differences, category-reference
     rebasing of scores and all curvature arrays, extreme/saturated tails,
     memory/numerical/input guards, deterministic RNG behavior, and fixed-design
     optimization/Hessian agreement with independent base-R references.
     No Gaussian/scalar-family fit route changed. This completes the internal
     observation-likelihood slice, not multinomial Laplace fitting.
   - **Conditional-mode validation:** isolated installation and the expanded
     focused regression set pass 1,974 assertions, zero failures/skips, with
     only the previously known structured-covariance false-convergence warning.
     Independent dense conditional optimization, numerical Hessians and
     determinant/likelihood calculations agree for coupled category designs
     and non-diagonal precision. Tests also cover scalar-binomial equivalence,
     multiple category/random-effect dimensions, category-reference rebasing,
     saturated curvature with symbolic reuse, evaluation-order invariance,
     iteration/overflow failure recovery, zero random designs, malformed or
     serialized contexts, invalid precision, workspace/index admission and RNG
     preservation. Base-only namespace loading and theory-article rendering
     also pass. Public nominal fitting and statistical calibration remain
     unqualified; the mode slice alone does not provide marginal block scores.
   - **Marginal-score validation:** isolated source installation and the
     expanded focused regression set pass 2,149 assertions, zero failures/skips,
     with only the known structured-covariance false-convergence warning.
     Independent dense implicit-mode/log-determinant derivatives and central
     differences of fully re-solved Laplace likelihoods agree for fixed effects,
     nonlinear Cholesky precision coordinates, and joint directions affecting
     both. Tests explicitly distinguish the full score from frozen-mode and
     omitted-determinant approximations. Binary/scalar equivalence, category
     rebasing, fixed-coordinate transformations, zero random designs,
     saturation, covariance proposal order/failure recovery, true sparsity
     changes and panels beyond 32 columns are covered. One-column-budget solves
     agree with default panels; disabled gradients perform no derivative solves.
     Analytic and numerical-score fixed-effect optimizers reach the same
     optimum. Input/workspace guards and RNG preservation pass, as do
     base-only namespace loading and theory-article rendering. These establish
     numerical derivatives, not nominal statistical coverage or a complete
     public fitted-object contract.
   - **Nominal preparation/coordinate validation:** isolated installation and
     focused nominal/multinomial/ordinal/Laplace/ZIP/bootstrap/PQL/Gaussian
     regressions plus existing modular-control/capability contracts pass
     2,508 assertions, zero failures/skips, with only the previously known
     structured-covariance false-convergence warning. Tests cover labelled
     reference coding, numeric/unsupported responses, common missing-data
     filtering and category retention, transformed/full-rank fixed bases,
     contrasts/dot expansion, offsets, shared/independent/correlated random
     designs and existing-control rejection. Independent Kronecker precision
     calculations with supplied nonidentity relationship precision match
     likelihoods, modes and Hessians. Fully re-solved finite differences agree
     with scores on native covariance coordinates for unstructured category
     covariance, multiple terms and fixed flags. Binary descriptor objectives
     agree with the scalar engine; zero-free-parameter cases, proposal/mode
     failure recovery, cache invalidation/result isolation, RNG preservation
     and joint fixed/covariance optimization pass. Base-only namespace loading
     and the beginner-friendly theory article render also pass. No shared
     installed library was replaced, and this validation does not qualify
     nominal coverage or a public fitted-object interface.

   - **Public nominal fit and downstream milestone:** exposed
     `multinomial(reference=)` for native Laplace ML and assembled a complete
     fitted object with labelled category probabilities, fixed/random
     estimates, likelihood, AIC/BIC and optimizer diagnostics. Conditional and
     zero-random-effect prediction support new data; probability uncertainty
     uses post-fit joint fixed/random curvature and a delta-method transform.
     `summary()`, `vcov()`, `postPEV()`, response/Pearson residuals, nested
     fixed-effect comparisons, factor-preserving conditional/population
     simulation, and nominal-specific parametric bootstrap refits are
     connected. Population-marginal prediction and deviance/working/quantile
     residuals remain unsupported and are rejected explicitly rather than
     routed through scalar-family code.
     Installation and focused nominal/multinomial/ordinal plus scalar-Laplace,
     bootstrap, PQL and prediction regressions pass; the established
     structured-covariance false-convergence warning remains. No package
     dependency or Gaussian/scalar fitting route changed. This does not
     qualify nominal statistical coverage.

The next deliverable is nominal statistical calibration and broader
simulation-diagnostic studies, followed by full frozen-snapshot package
qualification. Keep the current explicit boundary guards, post-fit
uncertainty contract, and restricted population-marginal covariance support.

For the bounded ZIP formula slice, the next statistical priority is low-count
identification/boundary diagnostics and calibrated inference: five registered
low-rate gates failed. Do not extend the zero component's model class first.

The release gate remains a stable, scientifically checked native likelihood
with explicit limitations, not just a successful optimizer return code.
