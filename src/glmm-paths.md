# Native GLMM feature implementation plan

Status: roadmap with workstream 2 implemented in the working tree, 2026-11-03.

This plan follows the current native GLMM implementation and the ordered
priorities agreed for the next feature work. It complements
[`native-glmm-epic.md`](./native-glmm-epic.md), which records implementation
and validation history. This document is a forward plan, not a claim that the
work below is implemented or scientifically qualified.

### Execution update (2026-11-03)

The first ZIP-diagnostics implementation slice is now in the working tree.
Native ZIP fits report per-trait observation/zero/positive-count totals,
zero-design rank, fitted structural-probability range, and expected
structural-versus-count-generated zero allocations among observed zeros.
Low-information flags are explicitly heuristic; diagnostics use accepted-fit
values only and do not compute Hessians, standard errors, or PEV during fit.
Focused ZIP-formula and native-family tests pass (419 assertions total, no
failures, warnings, or skips), and the GLMM guide renders with the diagnostic
description. This does not complete the ZIP workstream: diagnostic thresholds
and interpretation still need review, and sparse/mixed-family qualification
remains required.

Workstream 2 is implemented end-to-end for categorical `familym()` native
Laplace ML: ordinal+scalar, nominal+scalar and ordinal+nominal+scalar share an
observation map and native block likelihood. Numeric long-format category
codes retain labelled probability outputs. Fixed nominal logits expand by
category; random formulas do not silently create independent effects.
Explicit `.familyCoordinate` terms request separate/correlated coordinates;
ambiguous single-family `.nominalCategory` expansions error.

Mixed `fitted()`/response prediction return trait-named vectors/probability
matrices. Prediction (including conditional/full post-fit uncertainty),
simulation, residuals, summary/vcov, postPEV, nested fixed-effect comparisons
and boot_mmes are wired and tested. Nominal working/quantile residuals,
classify aggregation and population prediction SE remain genuinely
unsupported; population integration retains existing nominal intercept/
identity-relationship/category-count restrictions.

The scalar backend now admits full nominal observation-block curvature using
the existing softmax kernel. Scalar/ordinal outer scores remain analytic;
mixed nominal outer scores use central differences of the native joint
Laplace objective, including mode/determinant dependence. SE and covariance
work is post-fit only, with existing controls and memory policy. Gaussian
and PQL routes are unchanged. No candidate extension from workstream 3 was
implemented. The small R-independent `src/mme_glmm_paths.h` decision policy
is now wired into actual dispatch and execution-plan reporting, alongside
the likelihood-layout and gradient labels.

Validation: package native compilation/install, R/Rd syntax checks and both
GLMM article renders succeeded. Twelve focused scalar, categorical, ZIP,
PQL and bootstrap suites passed 2,128 assertions with zero failures; the
existing structured-random Laplace fixture emitted its optimizer
false-convergence warning. The new mixed-family suite had no warnings.

## Goals and invariants

1. Improve low-count ZIP identifiability and boundary diagnostics before
   widening the mixture model class.
2. Extend `familym()` so ordinal and nominal responses can be fitted alongside
   existing scalar families, with observation-level likelihood semantics.
3. Select one subsequent family extension from demonstrated user demand rather
   than implementing beta-binomial, zero-inflated NB2, or hurdle models by
   default.
4. Preserve distinct solver/model paths, existing Gaussian behavior, and
   explicit requested-versus-selected diagnostics.
5. Keep all standard errors, parameter covariance, and PEV work post-fit.

The existing ML/Laplace objective, conditional-mode solver, covariance
descriptor machinery, and `mmes_options`/`mmes_control()` contract remain the
foundation. Do not route GLMM derivatives through Gaussian AI updates or
modify Gaussian objective equations to implement these features.

## Current baseline and boundaries

- Scalar native Laplace families and mixed scalar `familym()` models are
  implemented. ZIP supports a fixed-effects zero formula but not zero-part
  random effects.
- Single-family ordinal and nominal models have public native Laplace routes.
  Both categorical forms now participate in numeric long-format `familym()`
  through an observation-aware native scalar/category-block response path.
- Laplace requires ML/Henderson, at least one random term, and independent
  observation likelihoods. Structured residual covariance and supplied
  observation weights remain outside the supported route.
- Laplace fitting defers parameter uncertainty and PEV. Use the existing
  post-fit `vcov()`, `summary()`, prediction uncertainty, and `postPEV()`
  surfaces where applicable.
- Existing controls include `task$glmm` for estimator and GLMM-specific
  settings, and `execution$denseMemoryMB` for workspace admission. New controls
  must be justified, validated, documented in the options catalog, and
  rejected on unsupported routes rather than silently ignored.

## Architecture: GLMM path planning

Implemented in `src/mme_glmm_paths.h`, the GLMM counterpart to
`src/mme_paths.h`. This is a pure, R-independent decision module: it selects
eligible likelihood/execution paths and stable reasons, without fitting,
evaluating likelihoods, choosing optimizer steps or computing uncertainty.
`src/mme_glmm_paths.cpp` exposes the policy through a thin Rcpp bridge.
`R/FUN_glmm_paths.R` maps existing family/control/formula capabilities into
the policy inputs; no new public controls are introduced.

Before calling a GLMM fitter, `R/mmes.R` checks the native plan and dispatches
using its **selected** path. Ineligible requests fail before either fitter
runs. Gaussian and internal preparation/working-model routes bypass this
GLMM-only policy. Existing backend guards remain defensive checks for direct
internal calls; the public decision is made by this one policy.

The policy distinguishes `pql`, `native-laplace-scalar` (including ordinal
scalar latent predictors) and `native-laplace-category-block` (any nominal
trait). Nominal likelihoods are never converted to independent scalar rows.
PQL preserves its existing scalar/structured-residual objective semantics;
categorical and ZIP PQL requests are rejected. Native Laplace requires ML
and independent observation likelihoods (`rcov` omitted or `~units`).

`executionPlan$glmm$path` stores the actual dispatch decision:
`requested`, logical `eligible`, `selected`, stable `reason`, explanatory
`message`, `likelihood` (`scalar` or `category-block`), categorical composition
and mixed-family status. `requestedPath`, `eligible`, `selectedPath` and
`reason` are also exposed directly under `executionPlan$glmm`; the existing
human-readable `selected="Native Laplace"`/`"PQL"` label is unchanged.
Rejections select `"none"` without fallback and carry the same decision in
`condition$glmmPathPlan` on a `sommer_glmm_path_error`.

Stable reason codes are `eligible`, `categorical-pql-unsupported`,
`mixture-pql-unsupported`, `laplace-requires-ml` and
`laplace-residual-dependence`. Other numerical, covariance, response-schema
and workspace checks stay in their owning modules.

Planner validation: native package compilation/install and R syntax checks
passed. Eight focused planner, mixed-family, scalar Laplace, PQL, ordinal,
nominal and control-contract suites passed 1,368 assertions with zero
failures. The existing structured-random Laplace fixture retained its
false-convergence warning. Tests cover actual mixed dispatch/reporting and
prove that rejected PQL/REML/structured-residual requests cannot invoke either
fitting backend.

### Responsibilities

- Represent the requested estimator, response-likelihood composition,
  conditional-mode backend, and requested output products separately.
- Check mathematical eligibility from capabilities (for example scalar
  observation likelihood versus category-block likelihood), not family names
  alone.
- Return requested, eligible, and selected paths with a stable rejection or
  fallback reason. Never silently change the objective or reinterpret a
  category-block likelihood as independent scalar rows.
- Keep workspace admission and numerical solver policy in their owning
  modules. The planner can consume an eligibility result or memory budget but
  should not duplicate sparse factorization, CHOLMOD lifecycle, or resource
  accounting.
- Make post-fit inference a distinct phase/product request. The fit plan must
  not trigger a Hessian, fixed/random PEV, or standard errors.

### Policy scope

Compact enums/records cover the current decision inputs and outputs:

- likelihood layout: scalar-per-observation or block-per-observation;
- family composition: scalar-only, scalar plus ordinal, scalar plus nominal,
  or scalar plus both categorical forms;
- mode and objective path: PQL, scalar native Laplace, or categorical-block
  native Laplace;
- status: eligible or unsupported, with no silent fallback;
- reason code: categorical/mixture PQL, ML objective or residual dependence.

This intentionally does not create a second generic solver framework or
move Gaussian `mme_paths.h` policy, response validation, or memory admission
into the GLMM header.

### Integration boundaries

- `R/FUN_mme_modules.R` remains responsible for task resolution and selected
  execution-plan reporting.
- `R/FUN_mme_control.R` and the generated `mmes_options` catalog remain the
  source of validated public controls and user-facing control descriptions.
- `R/mmes.R` owns top-level estimator dispatch using the native policy;
  `R/FUN_glmm.R` owns Laplace fitting and shared observation alignment.
- `R/FUN_glmm_families.R` owns family specifications and family-parameter
  mapping; ordinal and nominal category semantics remain in their specialized
  modules.
- `src/mme_glmm_family.h` and `src/mme_glmm_multinomial.h` remain likelihood
  kernels. Sparse mode, precision, and marginal score work stays in their
  current native modules or extracted GLMM-specific helpers.
- `R/FUN_glmm_ordinal.R` and `R/FUN_glmm_nominal.R` remain owners of
  family-specific coding, prediction, residuals, and fit metadata.

The capability matrix must be exercised by tests: every requested/selected
path and every unsupported combination needs an explicit result, rather than
testing only that a fit happens to succeed.

## Workstream 1: ZIP low-count identifiability and boundary diagnostics

### Objective

Make weakly identified ZIP fits recognizable and interpretable without
silently discarding fits, changing the likelihood, or claiming that optimizer
convergence establishes identifiability. Preserve the current constant and
fixed-effects zero-formula model class.

### Implementation steps

1. **Audit current diagnostics and objective metadata.** Inventory the
   zero-component design, fitted probabilities, count means, optimizer score,
   covariance boundaries, joint curvature gates, rejected proposals, and
   existing `laplaceDiagnostics`. Record which quantities are available at
   the final accepted fit without invoking uncertainty calculations.
2. **Add deterministic structural checks.** Report zero-design rank and
   constant/aliased columns, response zero fraction and positive-count
   information, extreme fitted zero probabilities, near-separated
   zero-component predictions, and the degree of overlap/confounding between
   count and zero components. Define checks using the retained ZIP-trait rows
   only in mixed scalar-family fits.
3. **Add likelihood-based weak-identification indicators.** Where cheap and
   numerically well-defined, record profile/curvature indicators based on the
   already available accepted-point objective and score. Do not form or invert
   the outer Hessian during fitting. Do not label a heuristic as a formal
   identifiability test.
4. **Separate statuses.** Report optimizer convergence, mode convergence,
   boundary state, zero-component separation/weak information, and
   post-fit-uncertainty availability as distinct fields. An unavailable
   standard error must remain unavailable; never substitute zero or infer
   reliability from convergence alone.
5. **Expose the findings consistently.** Add concise print/summary warnings
   and structured diagnostics on the fit. Preserve compatibility of existing
   `laplaceDiagnostics` fields. If a new user-selectable diagnostic policy is
   needed, put it under `task$glmm` with strict `mmes_control()` validation;
   prefer no new option if existing diagnostics can express the result.
6. **Document interpretation and next steps.** Explain low-count and
   near-boundary risks, recommended simpler alternatives, and the distinction
   between diagnostic flags and inferential guarantees.

### Acceptance and tests

- Exact/near constant zero designs, rank-deficient formulas, complete or
  near separation, low and moderate count rates, structural-zero probability
  near zero and near one, and near-zero random-intercept variance.
- Mixed `familym()` data where only one trait uses ZIP, including missing zero
  covariates and retained-row filtering.
- Deterministic diagnostic output across repeated starts and row order; no
  mutation of fit objects or random-number state.
- Diagnostics distinguish optimizer failure, a likelihood boundary, weak
  information, and a valid interior fit. Rank deficiency errors only where
  the model is algebraically invalid; weak-information cases remain explicit
  warnings/diagnostics rather than arbitrary new rejection rules.
- Mocked tests prove fitting does not call `optimHess()`, `vcov()`,
  `summary()` uncertainty, or `postPEV()` to create these diagnostics.
- Focused ZIP tests, mixed-family/PQL regressions, control/capability tests,
  documentation checks, then a frozen-source full package check before release.

### Exit gate

Do not begin a more complex ZIP-family expansion until diagnostics are
implemented, documented, and shown to behave on registered sparse and
moderate-information scenarios. A statistical screening study can assess
whether the diagnostics are useful, but a broader calibration claim remains
a separate qualification.

## Workstream 2: mixed ordinal and nominal `familym()` support

### Target behavior

Allow ordinal and nominal categorical likelihoods to participate in the
existing long-format `familym()` model alongside supported scalar families.
The first supported contract should include scalar-plus-ordinal and
scalar-plus-nominal; the unified architecture and final integration tests must
also cover ordinal and nominal traits coexisting with scalar traits in the
same fit. Support one or multiple declared categorical traits only where the
category-major design, indexing, memory, and identifiability contracts are
explicitly satisfied.

Likelihoods remain conditionally independent by observation given the shared
random effects. This work does not imply structured `rcov`, PQL support for
categorical families, or independent binary expansions of nominal responses.

### Implemented workstream contract

The steps below record the intended architecture. The native block mode and
public surfaces are implemented; analytic nominal mixed marginal scores
remain a performance follow-up, not a missing likelihood or method surface.
Current central-difference scores include the complete marginal objective.

1. **Family/index contract.**
   - Generalize `familym()` preparation from one scalar family index per row
     to retained observation blocks with a family descriptor and a category
     map where needed.
   - Keep response, trait, trials, offsets, missingness, and category levels
     aligned through one retained-row map. A nominal observation contributes
     exactly one joint categorical likelihood block; it must never be treated
     as `K-1` independent Bernoulli observations.
   - Validate that every retained row maps to exactly one declared trait and
     that category labels/levels and reference coding survive filtering.
2. **Likelihood interface.**
   - Define a common native family evaluation contract for scalar contributions
     and categorical block contributions: log likelihood, predictor score,
     observed curvature, required directional contractions, and structured
     failure reason.
   - Keep ordinal cutpoint coordinates and nominal reference-category
     coordinates family-specific. Preserve the full cross-category nominal
     curvature.
   - Ensure mixed marginal scores include conditional-mode and log-determinant
     dependence for fixed effects, covariance parameters, family parameters,
     and ordinal cutpoints.
3. **Design and random-effects assembly.**
   - Reuse the existing model setup and covariance descriptors. Make shared,
     trait-specific, and category-specific random effects explicit in formulas
     and metadata; do not silently rewrite user formulas.
   - Preserve relationship precision alignment and covariance-parameter
     constraints across stacked categorical designs. Validate trait-location
     identification for ordinal fits and reference constraints for nominal
     fits after the combined fixed design is assembled.
   - Establish memory admission before category-major expansion and sparse
     assembly using `execution$denseMemoryMB` and existing memory policy.
4. **Outer optimization and fit assembly.**
   - Optimize all supported fixed, covariance, family, and threshold
     coordinates against the same joint marginal likelihood and native first
     scores. Preserve the scalar-only path and numerical results.
   - Add combined fitted-object metadata and ensure `logLik`, AIC/BIC,
     `familyParameters`, coefficient tables, boundary diagnostics, and
     optimizer diagnostics count and label each trait's parameters correctly.
5. **Public methods.**
   - Add trait-aware `fitted()`, `predict()`, `simulate()`, `residuals()`,
     `summary()`, `vcov()`, `postPEV()`, nested fixed-effect comparisons, and
     `boot_mmes()` support, or reject a method combination explicitly until
     implemented.
   - Return labelled category probability matrices for categorical traits;
     do not return category codes as if they were scalar means.
   - Keep uncertainty and PEV calculations deferred. Fit-time output requests
     for inverse/PEV continue to follow the native Laplace post-fit contract.
6. **Control and path reporting.**
   - Extend the GLMM path planner to distinguish scalar family indexing from
     category-block likelihood evaluation and report why PQL or an unsupported
     route cannot be selected.
   - Reuse `task$glmm` and existing `familyParameters` for estimator and
     family-specific settings. Add controls only when a real choice is
     necessary; validate combinations centrally and list them in the generated
     `mmes_options` catalog.
   - Keep `executionPlan$glmm` consistent with actual selected likelihood,
     mode, and inference phases.

### Acceptance and tests

- Scalar-only `familym()` results remain unchanged within established
  tolerances; single-family ordinal and nominal public tests remain unchanged.
- Mixed scalar/ordinal and scalar/nominal fits, plus an ordinal+nominal+scalar
  fit, agree with independent likelihood/score calculations on small fixtures.
- Verify category permutation/reference rebasing, category labels, absent
  levels, missing responses/covariates, offsets, trials on other traits,
  `naMethod*` behavior, and non-prefix row omission.
- Verify full nominal cross-category curvature and joint score finite
  differences, ordinal cutpoint score transformations, and joint covariance
  parameter derivatives in mixed models.
- Test shared/correlated random effects, identity and supplied relationship
  precision, fixed covariance parameters, and unidentifiable trait/category
  designs. Unsupported structured residual covariance and PQL paths fail
  before expensive allocation.
- Memory admission occurs before large dense/categorical allocations and
  respects existing resource controls.
- Mocked tests prove fitting performs no Hessian/SE/PEV work. Post-fit
  uncertainty and PEV tests verify trait-specific parameter ordering and
  covariance dimensions.
- Round-trip serialization, prediction with known/unseen levels, seeded
  simulation/RNG restoration, residual type semantics, comparisons, bootstrap,
  print/summary, docs, namespace loading, and source-package checks.
- Add each mixed-family test selector to the focused native GLMM validation
  suite; run Gaussian and scalar GLMM regression suites to certify unchanged
  paths.

### Explicitly out of scope for this workstream

PQL with ordinal/nominal outcomes; structured residual dependence; automatic
category collapsing; arbitrary shared-coefficient or category-specific offset
syntax beyond the explicit design contract; new quadrature approximations;
and any inference guarantee not supported by separate qualification.

## Workstream 3: select one family extension using demand

Do not commit to beta-binomial, zero-inflated NB2, or hurdle likelihoods until
there is documented evidence of user need. Before implementation, collect
feature requests/use cases and compare:

| Candidate | User problem addressed | Main design/validation cost |
|---|---|---|
| Beta-binomial | Extra-binomial variation in grouped-binomial data | Trial-aware likelihood, concentration/overdispersion parameter, boundary behavior and interaction with mixed-family trials |
| Zero-inflated NB2 | Excess zeros plus overdispersed counts | Joint shape and mixture parameter identifiability, potentially severe low-count boundaries, ZIP diagnostic reuse |
| Hurdle | Separate zero occurrence and positive-count process | Two-part likelihood semantics, zero-truncated count normalization, prediction/simulation/residual contract |

### Selection gate

1. Record concrete user requests, data scale, current workaround, and expected
   frequency across users; distinguish anecdotal requests from a recurring
   need.
2. Assess likelihood identifiability and data requirements, not only
   implementation overlap. In particular, ZINB must inherit the ZIP
   low-count/boundary diagnostics before it is considered safe to fit.
3. Estimate the smallest complete public surface: constructor, `familym()`,
   parameters, prediction, simulation, residuals, comparison/bootstrap,
   documentation, tests, and post-fit uncertainty.
4. Select one candidate with a written rationale and acceptance contract.
   If demand is tied or insufficient, defer rather than implementing several
   partial families.

### Implementation contract after selection

- Add a native family descriptor and likelihood/score/curvature implementation
  without adding a fitting-package dependency.
- Keep nuisance parameters in explicit family coordinates and include their
  full marginal score dependence. Include parameter counts, natural-scale
  estimates, boundary diagnostics, and delta-method uncertainty only in
  post-fit methods.
- Reuse common response indexing, missing-data policy, prediction, simulation,
  and diagnostics helpers where semantics are identical; avoid family-specific
  duplicated fit pipelines.
- Gate optimizer starts and invalid-domain proposals explicitly. Do not hide
  expected numerical rejections as programming errors or return success-shaped
  defaults.
- Require independent density and derivative references, mixed-family
  regressions, boundary/sparsity tests, seeded simulation, user-facing
  documentation, and frozen-source package validation before public release.

## Controls, documentation, and quality policy

- Public estimation settings stay under `mmes_control(task=list(glmm=...))`;
  family-specific starting/fixed parameters stay in
  `task$glmm$familyParameters`. Shared resource settings stay under
  `execution`, particularly `denseMemoryMB`, memory policy, and BLAS threads.
- `mmes_options` remains a catalog of validated controls, not a second
  inference or eligibility mechanism. Every newly accepted field needs a
  constructor/validation rule, owner, default, effect, fallback/error behavior,
  example, and generated-catalog regression.
- Update the GLMM user/theory vignettes and direct Rd documentation in the
  same feature tranche. Keep `src/native-glmm-epic.md` as the implementation
  record; reconcile its status/date and roadmap only when its existing staged
  changes can be integrated without overwriting contributor work.
- Every stage must report requested/selected execution paths, boundary and
  convergence status, and unavailable inference explicitly.
- Keep all standard errors, parameter covariance, prediction uncertainty,
  and PEV strictly post-fit. Fit-time diagnostics may use accepted-point
  values, scores, and cheap structural checks, but must not calculate a
  Hessian or inverse as a side effect.
- No stage changes PQL defaults, Gaussian fit routes, Gaussian objective
  equations, or claims of covariance catalog support without dedicated tests.

## Recommended sequencing and release gates

1. Implement and document ZIP diagnostics; pass focused regressions and
   low-/moderate-information screening.
2. **Implemented:** pure GLMM path planner, native bridge, dispatch/reporting
   integration and supported/rejected capability tests.
3. Implement `familym()` category-block response indexing, optimizer assembly,
   downstream public methods, and post-fit inference for mixed ordinal/nominal
   fits; validate each slice while preserving scalar-only results.
4. Collect demand and select exactly one additional family. Implement and
   validate its full public contract.
5. Reconcile this plan with the implementation epic and current options
   catalog, run a clean frozen-source package check, then treat platform and
   scientific qualification as explicit release gates.
