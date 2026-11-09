# Nonlinear mixed-model effects: phased implementation plan

## Progress checkpoint (2026-11-09)

| Phase | Status | Current position |
|---|---|---|
| 0: Traditional curve contract | Complete for approved scope | Versioned exponential, logistic and saturation specifications; no layers or backpropagation. |
| 1: Nonlinear fixed mean, Gaussian ML | Implemented and validated | Exact Gaussian integration over additive random effects, analytic scores and serializable point prediction. |
| 2A/2B: One random curve parameter | Implemented and validated within stated limits | Laplace ML, exact conditional curvature, full analytic scores, optional independent additive effects and distinct conditional/zero/population means. |
| 2C: Covariance coverage and latent reduced rank | Implemented and validated for stated scope | Native latent `rrm()` maps for additive and curve effects; dedicated FA and representative covariance checks. |
| 3: All existing response families (expanded approval) | Implemented and validated for stated scope | Complete native family-aware Laplace scores, scalar/categorical/mixed fitting, link/response prediction, joint-prior population integration, diagnostic residuals and conditional/population simulation. |
| 4: Post-fit inference | Implemented and validated for stated scope | Observed parameter curvature, natural/group covariance, conditional-mode covariance versus empirical Bayes approximate PEV, prediction intervals and regular comparisons. All calculations remain explicitly post-fit. |
| 5: Scalability and release qualification | Implemented and locally validated; external release gates pending | Profiling, sparse-component admission, independent numerical references, scaling, full regressions and package/vignette checks pass locally. Cross-platform CI and PDF-toolchain qualification remain unexecuted. |
| 6: Multiple random curve parameters | Implemented and locally validated | Independent named formulas and a separate joint parameter-axis covariance block; existing vsm syntax, complete cross-parameter derivatives and post-fit-only uncertainty. |
| 7: Classical-model consolidation and release gate | Planned; next priority | Close portability gates, extend independent numerical/statistical qualification and freeze the supported classical contract before graph work. |
| 8: Declarative graph and derivative foundation | Proposed; blocked by Phase 7 | Modular smooth graph evaluator and native backpropagation prototype; no public network fitting until derivative and API gates pass. |
| 9: Fixed multi-layer predictors | Proposed; blocked by Phase 8 | Small smooth fixed networks plus ordinary additive random effects; Gaussian ML first. |
| 10: Random multi-layer parameters | Proposed; blocked by Phase 9 | Small declared Gaussian parameter deviations, exact conditional curvature and complete Laplace derivatives. |
| 11: Graph families, post-fit inference and release | Proposed; blocked by Phase 10 | Qualify native families, identifiable uncertainty products, all output surfaces and resource/portability gates. |

Phases 7--11 are roadmap updates, not implementation claims or authorization
to begin network development. The next implementation priority is Phase 7.
Layers and backpropagation remain unavailable until the classical stability
gate is met and the graph API/implementation scope receives separate approval.

Latest Phase 3 validation checkpoint: **5,163 regression assertions passed,
zero failures**, with one existing linear-Laplace optimizer warning at
`test-mmes-laplace.R:271` (false convergence (8)). The gate includes all
nonlinear phases, control/capability contracts, fixed-rank handling,
categorical/mixed and scalar native-family paths, PQL/direct PQL, structured
Gaussian residual regressions, shared latent processes and existing RR tests.
A clean native rebuild, source-package checks and examples passed with no
errors or warnings and two existing notes (installed size and compiled-code
inspection). Standalone labeled nominal and zero-inflated categorical-mixture
models passed conditional/link/population prediction and seeded population
simulation after serialization into a fresh R process. A SparseQR unsorted-R
block-copy assertion exposed by category-specific fixed effects was repaired
and covered by the fixed-rank and nonlinear gates.

Phases 3, 4 and 6 are complete within the boundaries below. **Phase 5 is
implemented and locally validated; external platform and PDF-toolchain release
qualification remain pending**.
No PEV or parameter-uncertainty products
are computed during nonlinear fitting, point prediction, generative simulation
or default summaries. Explicit uncertainty requests are separate post-fit
operations.

Dependency policy: neither `nlme` nor `lme4` may be a package dependency,
including `Suggests`, or a requirement of tests, benchmarks or CI. Earlier
reference-comparison results below are historical evidence only. Those calls
have been replaced by independent base-R dense Gaussian ML checks and
simulated-truth recovery checks; no sommer fitting path used either package.
The replacement qualification gate passed 64 assertions across seven cases,
without failures or errors, with a maximum case time of 1.948 seconds.
The dependency-free joint benchmark passed all three group sizes, with maximum
parameter recovery error 0.058842 and generating-mean prediction RMSE 0.044188
(thresholds 0.15 and 0.1). Historical suite totals below predate this replacement.

Phase 4 acceptance evidence:

- 283 dedicated inference assertions passed without failures or warnings:
  independent dense Gaussian likelihood curvature and conditional covariance,
  natural-scale transforms, mode-motion PEV, all native scalar links,
  estimated family tails, ordinal/nominal and mixed categorical shapes,
  zero-component regression, grouped trials, all three curve families,
  Gaussian and Poisson random curves, additive/curve blocks, constraints,
  known relationships and singular reduced-rank restoration.
- The affected nonlinear Gaussian/inference gate passed 1,529 assertions
  without failures or warnings. The coupled legacy gate passed 5,172 assertions
  without failures, with the same existing `test-mmes-laplace.R:271` optimizer
  warning. These gates overlap and are not additive.
- Sixteen simulated Gaussian datasets checked parameter-interval coverage;
  400 conditional response simulations checked future-response coverage.
  These validate small known cases, not blanket finite-sample calibration.
- Cached and cold-reconstructed nominal/zero-regression inference passed
  serialization into a fresh R process. Exact fit-time inference omission,
  derivative agreement, numerical-range, memory and unsupported-request guards
  were checked explicitly.
- Clean source-package installation, S3/usage/code checks and examples passed:
  zero errors, zero warnings and one installed-size note. Vignettes were
  excluded from rebuilding; changed documentation was updated consistently.

Phase 2C validation checkpoint: 2,568 regression assertions passed with
one existing linear-Laplace optimizer warning. Final installed-package nonlinear,
shared-latent and existing reduced-rank tests passed 1,789 assertions without
warnings. Package checks and examples passed with no errors or warnings and
two existing notes (installed size and compiled-code inspection). Both additive
and curve reduced-rank models with known relationships passed conditional and
population prediction checks after serialization into a fresh R process.
These results certify the tested scope, not every covariance constructor.

Current covariance coverage:

- Positive-definite native covariance descriptors with analytic derivatives
  reuse the existing machinery for additive effects and independent/joint
  random curve coordinates. Dedicated nonlinear checks cover identity, diagonal,
  unstructured, AR1, compound symmetry, heterogeneous/homogeneous FA and
  representative Kronecker products. Known relationships and constraints are
  tested. This is not a blanket claim for every covariance constructor.
- `fam()` with positive specific variances has dedicated likelihood, score,
  fitting, conditional/population prediction and serialization tests.
- Pure `rrm()` uses normalized native loadings and a proper Gaussian latent
  prior. Both additive and random-curve effects support rank-one/rank-two and
  representative multi-axis maps, without a ridge or substituted positive
  specific variance. Full-coordinate modes are restored; latent point modes
  and dimension metadata are retained separately.

The phase descriptions below retain broader roadmap acceptance criteria.
Completed status above applies only to the approved implemented scope.
Phase 6 extends the earlier single-coordinate boundary to independent/joint
curve parameters. Exact-zero overall random variance, unseen prediction groups,
structured observation residuals and additive/curve cross-covariance remain
unsupported.

### Phase 3 expanded implementation and acceptance gate

The approved scope includes every existing native response family and link:
binomial (logit/probit/cloglog), Poisson, NB2, Gamma (log), Gaussian (identity),
zero-inflated Poisson with an optional zero-component regression, ordinal
(logit/probit), nominal softmax, and rowwise `familym()` mixtures. The curve
is a **shared contribution to each link predictor**, including every nominal
nonreference logit. Category-/trait-specific linear effects and covariance
terms retain the existing formula conventions. This does not introduce
separate nonlinear parameter sets by category or trait.

1. Reuse native family likelihoods, normalized densities, block curvatures,
   family-parameter derivatives and categorical preparation. Preserve the
   exact Gaussian additive route. Otherwise use modular exact-Hessian
   conditional modes and complete analytic Laplace scores, including
   family parameters, cutpoint gaps, zero-design coefficients, latent RR
   loading maps, implicit mode motion and determinant derivatives.
2. Reuse `task$glmm$familyParameters` for fixed/estimated nuisance parameters;
   nonlinear outer/mode/quadrature choices stay in their existing nonlinear
   controls. Do not silently use PQL, numerical gradients or other GLMM
   fitting settings. Preserve observation masks through categorical stacking.
3. Implement scalar and category-aware link/response outputs, retained/new-data
   designs, serialized prediction and existing diagnostic residual conventions.
   Population prediction must integrate both additive predictors and random
   curve deviations against their joint prior; it is not inverse-link at zero
   modes. Bound quadrature workspace and fail if accuracy is not attained.
4. Implement conditional and population response simulation using existing
   family samplers, including reduced-rank and relationship covariance,
   without computing uncertainty or changing explicit-seed RNG state.
5. Validate every likelihood/link, conditional Hessians and full outer scores
   independently; add integration/limiting-case oracles, fitting and public
   method regressions, mixtures, missingness, grouped trials, zero regression,
   covariance constraints/RR, serialization and resource failure tests.
   Run previous nonlinear and linear-family regressions and package checks.

Acceptance status (update only after verification):

- [x] Complete native scores/Hessians for every family, blocks and nuisance parameters.
- [x] End-to-end fits and all point-output surfaces, including mixtures and RR.
- [x] True population integration and conditional/population simulation.
- [x] Serialization, missingness, numerical/reference checks and legacy regressions.
- [x] Documentation, capabilities and package checks.

Independent homogeneous residuals remain the nonlinear residual contract.
Structured Gaussian residuals, multiple random curve parameters, layers,
backpropagation and all post-fit uncertainty are outside this phase.

### Phase 2C implementation (complete)

1. Establish an explicit tested covariance support matrix. Add nonlinear
   objective, score and prediction checks for `fam()` and representative
   diagonal, unstructured, serial and Kronecker structures, with fixed/free
   parameters, supported constraints and known relationships.
2. Introduce a modular latent-coordinate map for `rrm()`: write the singular
   effect as `b = L v`, with a nonsingular Gaussian prior on `v`. Reuse the
   existing reduced-rank conventions and preserve scale/ordering constraints;
   never invert the singular full covariance or add an undocumented ridge.
3. Begin with additive `rrm()` effects alongside the traditional fixed curve.
   Then extend the same map to one random curve parameter. The curve must be
   differentiated through its latent loading map; retain exact conditional
   curvature and loading-dependent terms in the full Laplace score.
4. Preserve conditional, zero-deviation and prior-population prediction.
   Save the loading/relationship basis for serialization and new-data
   reconstruction. Population integration uses the implied marginal prior
   variance, not fitted-mode uncertainty.
5. Keep numerical choices in `mmes_control()`; add a control only when there
   is a genuine supported choice. Keep PEV, SE and all uncertainty post-fit.

Completed acceptance gate:

- [x] Tiny reduced-rank models match an independent latent-space likelihood,
  conditional mode/Hessian and complete outer-score reference.
- [x] Loading, covariance and residual scores pass finite differences, including
  the effect of parameter-dependent loadings on the mode and determinant.
- [x] Exact additive likelihood agrees with an independent full-coordinate
  singular-covariance calculation; a one-dimensional nonlinear latent case
  also agrees with numerical integration within its Laplace tolerance.
  no artificial variance floor changes the statistical model.
- [x] Prediction, serialization, relationships, constraints, rank failures and
  memory admission have explicit tests.
- [x] Existing Phase 1/2, Gaussian, PQL and linear-Laplace regressions still pass.

Phase 2C is complete for the tested scope. Phase 3 extends response families
with complete derivative checks rather than introducing layers or backpropagation.

## Approved initial scope

The initial implementation uses traditional parametric nonlinear fixed means:
exponential decay, logistic growth and saturation. Each has a native evaluator,
named parameter layout, explicit analytic coordinate derivatives and
serializable point prediction. Numerical choices belong to
`mmes_control(task=list(nonlinear=...))`. Layers, neural-network architectures
and backpropagation are deferred; they are not prerequisites for Phase 1.
The model is `y = offset + X beta + f(x; theta) + Z u + epsilon`.
Gaussian integration over the additive linear random effects is exact.
PEV and all uncertainty are strictly post-fit; nonlinear post-fit inference
is deferred rather than approximated with the existing linear formulas.

The broader graph/network sections below describe future extensions, not
the current public API or the first release acceptance gate.

## Phase 2 implemented scope

`nonlinear_curve(..., random=list(rate=~vsm(ism(subject), sigma2=0.04)))`
declares one Gaussian random curve coordinate. Any one existing curve parameter
is eligible: log deviations for positive parameters, identity for midpoint.
Optional ordinary additive terms remain independent covariance blocks.
One unit-incidence curve coordinate per observation is required; native shaped
covariances and positive-definite known relationship precision are reused.

Dedicated native derivative, mode and Laplace modules retain exact conditional
curvature, safeguarded deterministic multistart optimization, and complete
analytic Laplace scores including implicit mode motion. Phase 1 retains its
exact Gaussian route. Choices live in `task$nonlinear$mode` and
`task$nonlinear$population`, not model callbacks or hidden options.

Conditional, zero-deviation and prior-population point means are distinct.
Population integration uses one-dimensional Gaussian quadrature with checked
successive-order agreement, not mode covariance or parameter uncertainty.
New-data predictions require represented groups and a saved covariance basis.
Group natural parameters and coordinate modes are retained without pointers.
All uncertainty and PEV remain deferred post-fit.

Phase 2 boundaries (multiple curve parameters are subsequently added in Phase 6):
no correlated
additive/curve blocks, random slopes, warm starts, stochastic starts, or
arbitrary singular covariance or zero overall random variance. Built-in `rrm()`
is supported through the Phase 2C latent representation.
The existing `vsm()` positive-variance contract
is preserved; with additive effects retained, remove a curve-random declaration
to fit exact zero rather than silently replacing it with a small variance.
Fixed-only fits remain unsupported. Multistart is not a global-mode
certificate.

## Longer-term goal

Add a native, likelihood-based nonlinear predictor for mixed models. The
predictor must support learned nonlinear fixed effects and, in a later phase,
nonlinear random effects. After the classical stability gate, a small smooth
one-hidden-layer network such as

```
eta_i = tanh(x_i W + c) v + d + z_i u
```

is a future initial target, not a special-case implementation. Multiple
hidden layers must compose through the same graph contract. Other supported
activation functions and composable nonlinear mean functions should use the
same parameter, derivative, fitting, prediction, and validation machinery.

Here `W`, `c`, `v`, and `d` are fixed parameters; `u` denotes conventional
additive random effects. A random nonlinear component instead has a form such
as `eta_i = f(x_i; theta, u_group[i])`, with `u` assigned a declared Gaussian
covariance. Keep predictor nonlinearity distinct from a response-family link:
the existing GLMM families transform the linear predictor into response
means, while this feature makes the predictor itself nonlinear in estimated
parameters.

## Acceptance conditions

1. A user can fit a supported nonlinear fixed predictor with Gaussian
   observations and sommer random effects, and retrieve fitted parameters,
   likelihood, diagnostics, and predictions for training and new data.
2. A later supported random-nonlinear fit integrates its random parameters
   under their declared covariance; it does not treat them as fixed features.
3. Supported nonlinear graphs share one validated parameter layout and
   derivative contract. Unsupported functions, dimensions, and combinations
   fail with explicit errors.
4. Existing Gaussian, PQL, and native Laplace linear-predictor behavior
   remains unchanged unless the new model form is requested.
5. Tests compare reported objective and derivatives against independent
   small-model calculations, and verify recovery on simulated identifiable
   cases.

## Scope boundaries

- Start with Gaussian maximum likelihood. Conventional REML relies on a
  linear fixed-effects design and does not directly extend to arbitrary
  nonlinear fixed parameters. Do not label nonlinear fits REML.
- Treat Gaussian nonlinear mixed models and non-Gaussian GLMM response links
  as separate implementation steps. Reuse family likelihoods only after the
  nonlinear predictor derivatives are correct.
- Do not silently convert nonlinear random coefficients to fixed effects,
  linearize them once and call that joint estimation, or fall back to PQL.
- Do not add a mandatory neural-network dependency or execute arbitrary R
  callbacks inside the native optimizer.
- The first implementation need not provide every covariance structure,
  residual structure, optimizer, or inference option. Reject unsupported
  combinations explicitly and report selected paths.
- Do not promise unique network coefficients. Hidden-unit permutations,
  ReLU rescaling, and other symmetries can yield equivalent fitted values.

## Source review: pre-implementation fitting architecture

This historical review motivated Phases 0--6. Those phases now provide native
traditional-curve evaluators, derivatives and output methods; the table below
describes the original linear-engine seams, not a claim that nonlinear support
is still absent. Future graphs should reuse the implemented nonlinear modules
rather than build a second mixed-model fitting stack.

| Area | Current owner and constraint |
| --- | --- |
| Public fitting dispatch | `R/mmes.R`: `.mmes_fit()` selects Gaussian, PQL, or native Laplace paths. GLMM path selection and compatibility checks live in `R/FUN_glmm_paths.R` and `src/mme_glmm_paths.h`. |
| Formula preparation | `.mmes_fit()` builds sparse fixed design `X`, random design blocks `Z`, and covariance descriptors. Preparation-only output is reused by `R/FUN_glmm.R`. |
| Existing Gaussian solvers | `src/ai_mme_sp2.cpp`, `src/mme_direct.cpp`, and shared covariance/residual modules solve linear mixed-model equations. They assume fixed and random contributions are design matrices times coefficients. Do not distort these engines to fit the nonlinear mean. |
| Native likelihood | `R/FUN_glmm.R` prepares inputs, initializes fixed effects, and optimizes fixed and covariance parameters with `stats::nlminb`. `src/mme_glmm_laplace.cpp` owns random-mode evaluation, Laplace objective, and native score calculations. |
| Predictor assumption | `laplaceModeEvaluation()` constructs `eta = offset + X * fixedEffects + Z * mode`. Its mode Hessian and score use `Z` directly. `LaplaceObjective::gradient()` uses `X` directly and differentiates the existing linear-predictor Laplace objective. |
| Random covariance | `src/mme_glmm_precision.h` and existing covariance descriptors prepare random-effect precision, determinants, and derivatives. Preserve this machinery for Gaussian random nonlinear parameters where its Gaussian-prior assumptions apply. |
| Families and derivatives | `src/mme_glmm_family.h` evaluates response likelihoods, score, curvature, and higher derivatives with respect to the link-scale predictor. These are not derivatives with respect to nonlinear model parameters. |
| Controls and reports | `R/FUN_mme_modules.R` validates controls and builds `executionPlan`; `R/FUN_mme_capabilities.R` exposes control capabilities. New route choices must be validated and reported here, not hidden in the fitter. |
| Result methods and prediction | `R/FUN_methods_mmes.R`, `R/FUN_glmm_families.R`, `R/FUN_glmm_system.R`, and `R/predict_mmes.R` provide likelihood, summaries, simulation, uncertainty, and prediction. Laplace predictions rebuild fixed/random designs from saved terms and levels. Nonlinear models need their graph and parameter mapping retained and evaluated on `newdata`. |
| Existing regression coverage | `tests/testthat/test-mmes-laplace.R`, the family-specific Laplace tests, `tests/testthat/test-predict-mmes-ci.R`, and GLMM theory/article documentation cover the linear-predictor route and should remain regression gates. |
| Native registration | `src/mme_glmm_laplace.cpp` uses Rcpp exports; generated wrappers live in `R/RcppExports.R` and `src/RcppExports.cpp`. Follow repository generation/build practice when adding native entry points. |

The key design conclusion: a nonlinear predictor is not another activation in
`family=`. It requires a predictor evaluator and derivatives with respect to
estimated fixed and random parameters. It must own its likelihood calculations
without changing the meanings of existing `fixed`, `random`, or `family`
formulas.

## Proposed future graph model and API contract

Use one declarative predictor graph for fixed and random nonlinear components.
Do not settle public syntax until a small API prototype proves that formulas,
parameter names, dimensions, prediction, and serialization work together.
This proposal does not replace `nonlinear_curve()`, its independent/joint
specifications, or existing `fixed`, `random` and `vsm()` formula syntax.

The API should express, at minimum:

- Input columns or existing prepared terms.
- An ordered list of layers/functions, each with a declared output size and
  activation or mean function.
- Which parameter blocks are fixed and which are group-specific random
  parameters.
- Parameter names, initialization, and supported constraints.
- Optional conventional sommer random effects added to the nonlinear
  predictor.

The specification describes model structure only. Numerical choices must use
the existing `mmes_control()`/options and capability-validation approach.
Initialization, derivative strategy, optimization, convergence and memory
budgets must not become a second graph-specific numerical API. Public names
and any new control fields remain subject to the Phase 8 API gate.

Prefer a validated layer specification or a small formula DSL over arbitrary
closures. A closure-based API complicates native evaluation, differentiation,
reproducible optimization, prediction on new data, serialization, and
capability reporting. Keep the graph immutable during a fit; only its
parameter vector changes.

Use distinct notation and parameter names for weights, biases, and random
effects. The current conversation's `b` can mean a coefficient vector or
bias; API output must remove that ambiguity.

For group-specific network parameters, define whether `u_group` is a
deviation around shared fixed parameters or a complete group-specific
parameter. Prefer deviations, `theta_group = theta_fixed + u_group`, with
`u ~ N(0, G)`, because it defines a population mean and reuses sommer
covariance structures. Document any parameter transform before applying a
Gaussian prior: a Gaussian prior on unconstrained log-parameters does not
make the transformed natural-scale parameter Gaussian.

## Phases

### Phase 0 — model contract and numerical prototype

**Deliverable:** a bounded traditional curve specification and an agreed public
model equation/API before touching Gaussian solvers. For the approved initial
scope use `nonlinear_curve(input, curve)` with exponential, logistic and
saturation curves, natural-scale starts, log-positive parameter coordinates
and explicit analytic derivatives. Graph/layer prototypes are deferred.

1. Specify whether the nonlinear predictor returns the response mean or the
   link-scale predictor. For initial scope, return the link-scale predictor;
   the response family remains responsible for its own link.
2. Define stable parameter ordering, dimensions, names, bounds/transforms,
   initialization, missing-value handling, and serialization.
3. Define derivatives needed for outer optimization and conditional random
   modes. The evaluator contract should return predictor values, Jacobians
   with respect to fixed and random parameter blocks, and required curvature
   products without constructing observation-by-parameter dense arrays when
   avoidable.
4. Compare analytic or automatic derivatives against central differences on
   tiny inputs. Define ReLU's derivative convention at zero and test behavior
   immediately on both sides of the kink.
5. Establish identifiability warnings and deterministic initialization
   behavior. Do not imply that different optimizers or starts must return
   identical network coefficients.

**Gate:** approve equations, API shape, derivative contract, and supported
parameter constraints. No changes to existing model paths.

### Phase 1 — nonlinear fixed predictor with Gaussian ML

**Deliverable:** end-to-end fitting for a nonlinear fixed mean plus ordinary
linear Gaussian random effects, using the existing covariance descriptor
machinery. Allow a fixed-only nonlinear Gaussian model only if its likelihood,
dispersion, and result contract are explicit; it is not required for the first
mixed-model release.

1. Add a dedicated nonlinear model specification and preparation layer. Reuse
   existing observation filtering, offsets, random-effect preparation, and
   covariance descriptors where their contracts match.
2. Add native exponential-decay, logistic-growth and saturation evaluators
   behind one shared value/coordinate-derivative contract. Do not add layers
   or backpropagation in this slice.
3. Optimize nonlinear fixed parameters jointly with covariance and residual
   parameters against the Gaussian ML objective. Keep the existing random
   predictor linear in this phase, so its established conditional mode and
   Hessian remain valid.
4. Use a numerical outer gradient initially only if needed as a reference or
   fallback. It must be explicit, bounded by validated parameter/workspace
   limits, and reported. Do not claim it scales to large networks.
5. Store the curve specification, parameter layout and fitted parameters.
   Add new-data prediction by evaluating the same native curve, without
   manufacturing fixed-design columns for estimated nonlinear parameters.
6. Route the request through an explicit model-path planner. Unsupported
   controls must error before fitting; do not silently select a linear fit.

**Initial restrictions:** Gaussian ML, independent diagonal residuals,
   Henderson/native Laplace fitting, and conventional additive Gaussian
   random terms. Reject REML, PQL, direct fitting, structured residuals,
   unsupported family mixtures, and unsupported control combinations until
   each has a proved implementation.

**Tests and gate:**

- At an equivalent fixed nonlinear offset, match existing Gaussian ML
  likelihood and covariance estimates within numerical tolerance.
- Traditional curve values and full Gaussian scores match independent
  dense likelihood and finite-difference calculations on small datasets.
- Simulated nonlinear fixed effects recover predictions and parameters in an
  identifiable fixture.
- Finite-difference checks validate parameter scores; invalid shapes, NaNs,
  rank/identification issues, and optimizer failures produce explicit errors
  or diagnostics.
- Existing Gaussian, PQL, and linear Laplace tests remain unchanged and pass.
- `predict()` reproduces training fitted values and handles `newdata`;
  serialization and refitted-object prediction preserve the graph.

### Phase 2 — nonlinear random parameters

**Deliverable:** group-specific nonlinear parameters integrated as Gaussian
random effects, not treated as a fixed design matrix.

1. Extend the conditional mode evaluator to accept `eta = f(X; theta, u)` and
   calculate the random score and Hessian. For nonlinear `f`, the exact
   Hessian includes both the likelihood-curvature/Jacobian product and
   second derivatives of `eta` with respect to `u`; `Z' W Z + G^-1` alone is
   only valid for a linear predictor.
2. Compute the Laplace determinant from the actual conditional negative
   log-posterior Hessian. Reject non-positive-definite or invalid curvature
   with a clear numerical reason.
3. Reuse covariance precision, determinant, constraints, and relationship
   preparation only where the random parameter ordering maps correctly to
   existing covariance descriptors. Define the mapping from network
   parameter blocks to groups and covariance factors explicitly.
4. Include all learned nonlinear fixed parameters, random-mode coordinates,
   and covariance parameters in the likelihood and derivative bookkeeping.
   Validate the outer gradient against finite differences of the complete
   Laplace likelihood, including mode and log-determinant dependence.
5. Define mode starts, warm-start validity, boundary diagnostics, and
   convergence reporting for non-concave modes. Warm starts must not conceal
   failed cold starts or change the selected objective.
6. Restrict initial random nonlinear effects to small group-level parameter
   blocks. Per-observation random network weights are an explicit non-goal
   until sparse structure and memory behavior are demonstrated.

**Tests and gate:**

- A Gaussian random-intercept/slope special case reproduces the existing
  linear Laplace likelihood and modes.
- A simulated nonlinear random-intercept or one-parameter random curve
  recovers population and group-level predictions.
- Analytic/automatic mode score and Hessian agree with finite differences.
- The full Laplace score agrees with finite differences at interior parameter
  points and across covariance parameterizations.
- Tests cover multiple groups, one group, zero variance boundaries,
  separation-like data, ReLU knots, invalid curvature, and deterministic
  repeated fitting.

### Phase 3 — additional response families

**Deliverable:** nonlinear predictors combined with already-supported
observation families when likelihood and derivative checks pass.

1. Connect the predictor evaluator output to `GLMMFamilyModel`'s response
   likelihood evaluation. Keep response-link derivatives separate from
   predictor-parameter derivatives.
2. Add family support one family at a time, beginning with a scalar family
   whose likelihood and derivatives have independent test oracles.
3. Re-derive outer scores for nonlinear fixed and random parameters. Existing
   `LaplaceObjective::gradient()` uses linear design matrices and is not
   directly reusable for nonlinear parameter blocks.
4. The expanded approval includes categorical block likelihoods,
   zero-inflated models and family mixtures. Validate their cross-derivative,
   indexing, prediction and simulation contracts in this phase rather than
   deferring them. Structured Gaussian observation residual blocks remain out;
   they are distinct from supporting Gaussian observations within a mixture.

**Gate:** each enabled family must pass likelihood, score, Hessian/mode,
prediction, simulation, and diagnostics tests. Unsupported family/model
combinations must fail explicitly.

### Phase 4 — inference and full sommer surface

**Deliverable:** defined uncertainty and user-facing behavior for nonlinear
model parameters and fitted responses.

1. Add post-fit observed-curvature or clearly labeled approximate parameter
   covariance. Transform uncertainty to natural-scale parameters with a
   documented delta method where appropriate.
2. Add prediction uncertainty using the nonlinear predictor Jacobian and
   supported random-effect covariance. Distinguish conditional, zero-random-
   effect, and population-marginal predictions; nonlinear transformations
   mean these are generally different.
3. Add `summary()`, `coef()`, `vcov()`, `fitted()`, `residuals()`,
   `simulate()`, `logLik()`, AIC/BIC, and model comparison only where their
   parameter counts and likelihood semantics are well-defined.
4. Add selection/inference for nonlinear terms only with a valid null
   hypothesis and boundary-aware calibration. Do not advertise ordinary Wald
   or likelihood-ratio tests for weakly identified or boundary parameters
   without evidence.
5. Update `R/FUN_mme_capabilities.R`, control documentation, `man/mmes.Rd`,
   and the GLMM/model-theory vignettes. Explain parameter symmetries,
   non-convex optimization, activation kinks, initialization, and limits.

**Gate:** all exposed methods either return mathematically defined outputs or
reject unsupported requests explicitly. Validate prediction and uncertainty
against simulation on small known models.

Traditional-curve implementation:

- `post_nonlinear()` reconstructs operation-local, serializable exact Gaussian
  or Laplace likelihood inputs; checked differences of analytic scores provide
  observed marginal parameter curvature. Native conditional curvature remains
  exact. No fitting-time Hessian or inverse is introduced.
- `task$nonlinear$inference` in `mmes_control()` owns relative difference step,
  successive-step agreement, standardized eigenvalue tolerance and conditioning
  thresholds. Memory admission uses `execution$denseMemoryMB`.
- `vcov()` distinguishes natural mean/all parameters, working coordinates,
  conditional-mode covariance and approximate random PEV. The empirical Bayes
  PEV adds estimated-parameter mode sensitivity and is explicitly not claimed
  to be an exact frequentist PEV. Constraints, nuisance parameters, moving RR
  maps and group-parameter transformations are propagated.
- `postPEV()` dispatches to this implementation. Mode 1 retains diagonal
  approximate PEV; mode 2 retains full random and parameter/mode covariance.
  Both currently reconstruct full conditional covariance under memory admission.
  No incompatible linear `Ci` is created. `postPEV(mode=0)` clears saved
  nonlinear uncertainty without computation.
- `confint()` transforms working-normal intervals for positive curve
  parameters and exposes working/all-parameter intervals. Descriptor reporting
  SEs reuse existing callbacks; `covmatrix_mmes()` and `vcDelta()` propagate
  the same covariance to matrix entries and scalar parameter transforms. Default
  summaries display saved inference only; `inference=TRUE` explicitly requests it.
- Prediction preserves scalar/matrix/trait-list shapes. Conditional covariance
  propagates both conditional random curvature and total parameter derivatives;
  zero-deviation and prior-population mean uncertainty propagate marginal
  parameters only. Scalar future-response normal intervals additionally include
  observation variability and, for population prediction, joint-prior variance.
  Discrete intervals are labeled approximations. Numeric future-category
  intervals are rejected; categorical probability confidence intervals are
  supported.
- AIC/BIC count all free optimizer coordinates. Comparisons are restricted to
  regular nested ordinary fixed terms with unchanged nonlinear architecture,
  data and random/covariance semantics. Explicit finite working-coordinate
  nulls permit approximate fixed/curve Wald tests. Automatic nonlinear-term
  selection, variance-boundary tests and unmodified linear-model uncertainty
  or reliability formulas are explicitly unsupported, not silently calibrated.
- Neural activations, kinks, layer symmetries and backpropagation remain outside
  the approved traditional-curve scope. Current diagnostics instead address
  local identifiability, nonconvex starts, convergence and curvature.

Completed acceptance gate:

- [x] Independent likelihood-curvature and conditional-mode covariance oracles.
- [x] Parameter, group and predictor uncertainty with explicit estimands/scales.
- [x] All existing native family/link and categorical/mixed output shapes.
- [x] Constraints, known relationships and latent reduced-rank restoration.
- [x] Regular comparisons and explicit rejection of uncalibrated boundary tests.
- [x] Simulation checks, serialization, resource guards and coupled regressions.
- [x] Documentation, capabilities, clean package checks and examples.

### Phase 5 — scalability and release qualification

This approved-scope step qualifies the traditional curves and
native families implemented above. Hidden-unit, graph and layer dimensions in
the broader roadmap below remain deferred until separately approved.

**Deliverable:** a qualified feature envelope, not a blanket claim of general
neural-network support.

Implementation and local qualification:

- Opt-in fitting timings use `diagnostics$profile="force"`; the separate
  `profile_nonlinear()` operation reconstructs native contexts without fitting
  or computing uncertainty unless explicitly requested.
- Sparse-support component admission precedes Hessian assembly and
  factorization, including categorical curvature and random-curve incidence.
  It is an operation-local conservative bound, not a process-wide RSS cap.
  Independent 10,000-coordinate admission and connected-fill rejection are
  checked directly without constructing a dense matrix.
- The seven new automated cases passed 62 assertions without failures,
  errors or warnings. Measured case times were 0.073--2.015 seconds on the
  local configuration, below the repository's 30-second per-case limit.
  At that checkpoint, independent `nlme` references covered exact additive
  Gaussian ML for all three curves and random-rate recovery/prediction.
  These have since been replaced by dependency-free dense ML and simulated
  truth checks, as recorded in the dependency-policy checkpoint above.
- The standalone benchmark passed all 18 cases: exact Gaussian additive,
  Gaussian random-rate and Poisson random-rate paths; 128--8,192 observations
  and 16--512 groups. All fits converged, recovered the generating curve
  within the declared absolute 0.35 tolerance (maximum error 0.2261), kept
  sparse designs and omitted fit-time uncertainty. Explicit RNG state was
  restored.
- At 8,192 observations and 512 groups, whole-fit elapsed times were
  4.582, 8.855 and 3.358 seconds respectively. Retained sparse design size
  was 660,992 bytes, and retained fitted-object sizes were
  1,627,552--2,306,400 bytes. These are inclusive descriptive timings and R
  object sizes, not disjoint phase costs or measured native peak RSS.
- Configuration: Linux x86_64, R 4.2.2, `g++ -std=gnu++14`, FlexiBLAS 3.2.1,
  one requested BLAS thread. Timings are not cross-platform guarantees.
  Generated Rcpp wrappers reproduced deterministically.
- Full source-tree regression coverage, followed by focused runtime-policy
  and non-CRAN Wald reruns, passed 16,351 checks across 743 cases with no
  failures or errors. There is one existing linear-Laplace optimizer warning
  (`test-mmes-laplace.R:271`) and one GPU-only skip on this CPU environment.
  The totals replace overlapping reruns rather than double-counting them.
- Two earlier nonlinear inference cases exceeded 30 seconds when measured
  across multiple families. Their independent family checks are now separate
  cases with unchanged numerical coverage. The complete nonlinear gate passed
  2,296 assertions across 89 cases without failures or warnings; every case
  was below 30 seconds (maximum 23.371 seconds). Three pre-existing,
  non-nonlinear cases exceeded 30 seconds in the full run; they are outside
  these additions and were not changed.
- A clean source-package installation, code/S3/Rd checks, examples, all
  vignette builds, vignette-code execution and vignette rebuilding passed:
  no errors or warnings and one installed-size note (16.9 MB). Tests were
  run independently, not skipped as a substitute for regression coverage.
- The CI matrix includes Linux release/oldrel/devel, Windows and macOS;
  those external jobs have not been executed in this session. PDF manual
  qualification remains pending because the local TeX/PDF tools are absent.
  Neither the historical minimum R version nor every compiler/dependency
  combination is certified by this local run.

Acceptance checkpoint:

- [x] Opt-in modular fitting and post-fit profiling, without implicit PEV.
- [x] Sparse-component admission and exact memory-threshold rejection tests.
- [x] Independent Gaussian references and simulated-truth scaling workloads.
- [x] Per-case nonlinear runtime policy, full regressions and deterministic wrappers.
- [x] Clean local package installation/check, examples and generated vignettes.
- [ ] Execute the external R/platform CI matrix before portable release.
- [ ] Complete PDF manual/toolchain release checks in a suitable environment.

1. Profile forward evaluation, derivative products, mode iterations, outer
   optimization, and memory separately across observations, groups, hidden
   units, and random parameter count.
2. Keep sparse mixed-model algebra for conventional random effects; avoid
   densifying `X`, `Z`, or covariance precision just to evaluate nonlinear
   terms. Bound any dense network work with existing memory controls or a
   model-specific validated budget.
3. Add scaling thresholds and explicit errors before allocations that exceed
   them. Report actual path, gradient method, iteration counts, starts, and
   convergence/boundary status in `executionPlan` and fit diagnostics.
4. Compare supported cases with an independent implementation and simulated
   truth. Record when estimates can differ because of initialization,
   optimizer, parameter symmetries, or Laplace approximation.
5. Qualify supported R versions, platforms, compiler behavior, package
   install/check, generated Rcpp wrappers, documentation build, and the full
   regression suite before release.
6. Keep each automated test case below 30 seconds and measure the new cases
   individually. Run slower recovery and scaling examples separately through
   `inst/benchmarks/nonlinear-qualification.R`, not as automated tests.

**Stop condition:** stop expanding when the documented supported envelope
passes numerical, regression, prediction, memory, and package qualification
gates. Add new activations or model structures only with their own derivative
and recovery tests.

### Phase 6 — multiple random curve parameters

Approved scope: all parameters of the existing exponential, logistic and
saturation curves may have Gaussian coordinate deviations, independently or
jointly correlated. Positive parameters retain log-coordinate deviations;
midpoint retains identity-coordinate deviations. No layers or backpropagation.

The existing syntax remains unchanged:

```r
nonlinear_curve("x", random=list(
  asymptote=~vsm(ism(group)),
  rate=~vsm(ism(group))))
```

Correlated parameters use a separate joint specification:

```r
nonlinear_curve("x", joint=list(
  parameters=c("asymptote","rate"),
  random=~vsm(usm(.parameter),ism(group))))
```

`.parameter` is an internally supplied covariance axis in the declared
parameter order, not a user data column or a change to `vsm()`. Existing
covariance descriptors, known relationships and constraints retain their
meaning. Parameters cannot occur in both independent and joint blocks.
One joint block may coexist with independent deviations on other parameters
and ordinary additive random terms.
The `.parameter` axis must be one covariance-shaping factor, never the final
grouping factor or a data column. Responses are not replicated or stacked.
Covariance reports qualify otherwise identical terms by `curve[parameter]`
or the joint parameter names; existing constraint indices retain their meaning.

Implementation and acceptance:

1. Version and validate the specification while retaining existing serialized
   single-parameter models and formula environments.
2. Carry explicit full-coordinate curve maps through preparation and moving
   reduced-rank latent designs; retain sparse mixed-model algebra.
3. Include every cross-parameter second/third derivative in exact conditional
   curvature and complete analytic Laplace gradients for all native families.
4. Wire joint conditional/zero/population prediction, simulation, profiling
   and cold-reconstructed serialized contexts.
5. Propagate joint group-parameter covariance, conditional curvature,
   empirical Bayes approximate PEV and prediction uncertainty only post-fit.
6. Validate independent dense derivative/likelihood references, known
   covariance and reduced-rank restoration, family/categorical shapes and
   simulation recovery. Keep each automated case below 30 seconds; move
   larger qualification workloads to standalone benchmarks.
7. Update public documentation, capabilities and examples and run coupled
   regressions and clean package checks before marking local completion.

Implemented architecture:

- Version-3 specifications retain old version-1/version-2 contracts. Explicit
  full-coordinate membership and group maps distinguish independent blocks,
  correlated parameter axes and ordinary additive effects.
- Sparse parameter designs are masked before the current reduced-rank map is
  applied. Analytic scores include moving loadings, implicit mode motion and
  every cross-parameter second/third derivative, including categorical
  likelihood blocks and family nuisance coordinates.
- Multiple-coordinate Gaussian fits reuse the exact-normal family-aware
  Laplace kernel; existing exact additive and optimized single-coordinate
  routes are unchanged. This is shared likelihood machinery, not a fallback
  estimator or Gauss-Newton approximation.
- Conditional, zero and joint-prior population prediction, simulation,
  profiling and serialized reconstruction all use the same maps. Positive
  coordinates retain log deviations and midpoint retains identity deviations.
- `vcov(fit,type="group")` explicitly requests the full natural-scale
  group-parameter covariance, ordered like `as.numeric(groupParameters)`.
  It combines parameter uncertainty, mode motion and conditional latent
  covariance. All SE, PEV, covariance and interval work remains post-fit.

Local acceptance evidence:

- The dedicated final joint gate passed 294 assertions across 31 cases,
  without failures, errors or warnings; maximum measured case time was
  10.402 seconds. Independent dense likelihood/Hessian references cover all
  three curves. Full numerical-score references cover every native family
  and link, mixed categorical blocks, zero regression and RR loadings.
- The final complete nonlinear/control gate passed 2,757 assertions across
  129 cases without failures, errors, warnings or skips. Every measured case
  was below 30 seconds; maximum time was 23.828 seconds. This gate includes
  the dedicated joint cases above; overlapping totals are not additive.
- Acceptance also covers independent-versus-diagonal-joint equivalence,
  reversed parameter order, all three logistic coordinates, joint plus
  independent/additive effects, known relationships, constraints, full/latent
  restoration, distinct covariance names and an independent post-fit
  group-covariance oracle.
- Coupled native-family, categorical, shared-latent, fixed-rank, reduced-rank,
  backend, acceleration and Wald regressions passed 2,548 assertions with no
  failures or errors. One existing linear-Laplace optimizer warning remains
  (`test-mmes-laplace.R:271`, false convergence (8)).
- A fresh R process passed cold/cached nominal joint inference, conditional
  and population prediction, intervals and seeded population simulation.
- The standalone `inst/benchmarks/nonlinear-joint-qualification.R` passed
  correlated asymptote/rate recovery at 8, 16 and 32 groups (16 points each).
  Maximum natural-parameter recovery error was 0.058842 against the declared
  0.15 tolerance. Maximum parameter and conditional-prediction differences
  against historical independent `nlme::nlme` fits were 0.003383 and 0.002205, below
  declared 0.08 and 0.1 tolerances. The approximations differ: likelihood
  equality is not asserted. RNG state was restored, designs stayed sparse,
  and no fit-time uncertainty was computed.
- Clean final-source package installation, code/S3/Rd checks, examples,
  all 13 vignette builds, vignette-code execution and vignette rebuilding
  passed: zero errors, zero warnings and one installed-size note (19.4 MB).
  Tests ran independently through the gates above, not as part of this
  `--no-tests --no-manual` check. The package-building tar utility separately
  normalized an out-of-range account UID to `nobody`.

Acceptance checkpoint:

- [x] Compatible independent/joint syntax and explicit full-coordinate maps.
- [x] Complete cross-parameter likelihood derivatives and reduced-rank maps.
- [x] All point-output and post-fit uncertainty surfaces, including serialization.
- [x] Independent numerical references, recovery and coupled legacy regressions.
- [x] Final complete nonlinear/control per-case runtime gate.
- [x] Clean local package installation/check, examples and all vignettes.

The external platform CI matrix and PDF-toolchain release gates listed in
Phase 5 remain pending; local completion does not certify those environments.

### Phase 7 — classical-model consolidation and stability gate

**Priority:** complete this phase before implementing layers or
backpropagation. "Rock solid" means a reproducible, explicitly bounded
classical support envelope, not every possible model combination.

**Deliverable:** a release-qualified classical baseline and an evidence-backed
decision on whether graph development may begin.

1. Close the outstanding Phase 5 external R/platform CI and PDF/manual gates
   against the Phase 6 code. Re-run the full repository suite, not only
   nonlinear/coupled subsets, and distinguish existing warnings/skips from
   regressions. Record versions, reproducible commands and evidence.
2. Maintain a support matrix crossing curve, family/link, independent/joint
   parameter layout, covariance descriptor, known relationships, constraints
   and reduced-rank maps. Use representative interactions and risk-based
   coverage; mark untested combinations rather than imply an exhaustive
   Cartesian-product guarantee.
3. Stress complete likelihood scores, exact mode Hessians and post-fit
   Jacobians with independent small-model oracles at multiple interior points.
   Include unequal group sizes, missing/filtering patterns, extreme but valid
   inputs, parameter-order changes, covariance conditioning and moving RR maps.
   Check objective/score consistency at the actual accepted optimizer point.
4. Qualify statistical behavior in standalone repeated simulations: bias and
   recovery, convergence/mode-selection rates, parameter/prediction interval
   coverage, covariance estimation and population quadrature accuracy.
   Predeclare tolerances and Monte Carlo uncertainty. Compare independent
   implementations where estimands/approximations are comparable; do not
   require equal likelihoods for different approximations.
5. Exercise weak identification, nearly zero variance, poorly conditioned
   covariance, non-positive final curvature, integration nonconvergence and
   memory rejection. Each must give an actionable diagnostic or explicit
   unsupported error, never a silent estimator switch, ridge or fabricated SE.
   Audit deterministic starts and repeated fits without claiming global modes.
6. Measure fit and explicit post-fit costs separately over observations,
   groups and random-coordinate count. Include joint/RR/categorical workloads,
   conservative workspace admission and standalone peak-memory measurement
   where available. Do not equate retained R object sizes with peak native RSS.
7. Freeze regression fixtures for old serialized specifications, formulas,
   named constraints and every public output shape. Test cold reconstruction
   in fresh processes and ensure default fits/summaries/predictions still do
   not compute PEV or parameter uncertainty.
8. Reconcile documentation, capabilities and examples with the tested matrix.
   Record remaining limitations and any migration policy before declaring
   classical stability. Keep each automated case below 30 seconds; use
   standalone scripts for longer recovery, coverage and scaling experiments.

**Exit gate:**

- [ ] Full current-source regressions and external platform/manual gates pass,
  or any platform exclusion is explicitly justified and approved.
- [ ] Numerical, statistical and resource checks meet predeclared thresholds;
  no unresolved task-caused correctness defects remain in the supported scope.
- [ ] Public compatibility, serialization and post-fit-only inference hold.
- [ ] Supported envelope and residual risks are documented.
- [ ] User approves proceeding to the graph/API prototype.

#### Classical extensions backlog — separate scope decisions

These are real gaps, not automatically part of the stability gate or promises
of current support. Prioritize them by user need after consolidation; each
requires its own design, implementation and qualification before admission:

- Fixed-only nonlinear models, without a dummy random term.
- New/unseen-group prediction and simulation, with explicit conditional versus
  prior-population estimands and relationship-based conditioning where valid.
- Exact-zero variance and boundary-aware model comparisons; never approximate
  a requested zero by a small positive variance.
- Structured Gaussian observation residuals, distinct from link-scale random
  dependence, with valid determinants, derivatives and prediction semantics.
- Cross-covariance between additive and curve parameters, non-unit random
  loadings/slopes, or more than one joint block.
- Trait/category-specific curve parameter sets rather than the current shared
  link-predictor contribution.
- Larger post-fit problems via selected solves/contrasts instead of full
  covariance reconstruction; preserve mathematical equivalence and memory
  admission, and do not move uncertainty into the fitting loop.
- Additional classical curves and boundary/weak-identification inference only
  with independent derivative, recovery and calibration evidence.

None of these must be bundled into network work. Existing rejections stay in
place until their respective support gates pass.

### Phase 8 — declarative graph and backpropagation foundation

**Prerequisites:** Phase 7 exit gate and explicit approval of this scope.
**Deliverable:** a bounded native graph/derivative prototype, not a new
production estimator.

1. Prototype a serializable, versioned graph with named input transforms,
   affine blocks, smooth activations, dimensions and fixed/random parameter
   maps. Start with identity and `tanh`; admit other smooth operations only
   with stable values and derivatives through the required order. Reject
   arbitrary callbacks, cycles and unsupported operations.
2. Keep graph validation/layout, forward evaluation, differentiation,
   likelihood/mode solving and post-fit output in separate modules. Reuse
   existing covariance descriptors, family kernels and factor/solve machinery.
   Preserve the traditional-curve API and optimized routes.
3. Implement native reverse-mode differentiation (backpropagation) for
   first derivatives of the declared graph. Backpropagation is a derivative
   calculation, not an optimizer, a response link or a replacement for
   mixed-model likelihood integration. Classical curves do not require layers.
4. Specify and validate Jacobian products, Hessian-vector products and the
   higher-order directional contractions required by exact nonlinear random
   curvature and complete Laplace scores. First-order backpropagation alone
   is insufficient for random network parameters. Avoid allocating dense
   third-order parameter tensors.
5. Check native derivatives against independent analytic tiny graphs and
   finite differences, including cross-layer paths, shared parameter indexing
   and transformed coordinates. An identity graph must reduce to the linear
   predictor; supported curve compositions must agree with the classical
   evaluator without rerouting existing curve fits.
6. Prototype new-data evaluation and serialization before settling public
   syntax. Store any input scaling learned from training data and reuse it
   unchanged. Numerical choices stay in validated `mmes_control()` options;
   runtime graph/tape memory is subject to existing execution budgets.

**Exit gate:** approved graph API; exact dimensions/ordering; independently
verified first/higher derivative products; serialization/new-data parity;
bounded memory; unchanged classical regression suite. No network fitting is
advertised until Phase 9 passes.

### Phase 9 — smooth fixed multi-layer predictors

**Prerequisite:** Phase 8. **Deliverable:** opt-in Gaussian ML with a small
smooth fixed network and conventional additive Gaussian random effects.

1. Begin with one hidden layer, then require a two-hidden-layer acceptance
   fixture to prove true composition rather than a hard-coded network shape.
   Keep initial widths/depth bounded by demonstrated resource behavior.
2. Reuse exact Gaussian integration over additive random effects. Optimize
   network and covariance parameters jointly with verified likelihood scores;
   do not replace ML with a squared-error training loss or silently add PQL.
   Keep the established `nlminb` route unless a separate optimizer is approved.
3. Define reproducible scaled starts, optional seeded initialization and any
   explicit multistart policy through the controls contract. Report chosen
   start, objective, convergence, gradient and competing-start outcomes.
   Preserve the caller's RNG state and never promise a global optimum.
4. Address hidden-unit permutations, scaling symmetries and confounding with
   ordinary fixed terms. Reject structurally redundant layouts where feasible;
   otherwise diagnose nonidentification. Recovery tests compare predictions
   or identifiable functionals, not arbitrary equivalent weight vectors.
5. Wire preparation, coefficient names, likelihood, diagnostics, prediction,
   residuals, simulation and profiling. Save all graph/preprocessing metadata,
   never operation-local pointers. No fit-time PEV or parameter covariance.

**Exit gate:** linear/curve special-case parity, independent Gaussian objective
and score references, identifiable simulated recovery, deterministic repeats,
training/new-data and serialization agreement, memory guards and full classical
regressions. ReLU, dropout, batch normalization, stochastic/minibatch training
and GPU execution are not implied by this phase.

### Phase 10 — small random multi-layer parameter blocks

**Prerequisite:** Phase 9. **Deliverable:** explicitly declared Gaussian
deviations on selected smooth network parameters, starting with small
group-level blocks rather than every weight or per-observation networks.

1. Preserve existing `vsm()` formulas and declare parameter-to-covariance maps
   separately, following the classical independent/joint contract. Establish
   fixed-plus-deviation coordinates, transforms and identifiability before
   admitting a block; reuse constraints, relationships and native RR maps.
2. Form the exact conditional random-mode Hessian, including score-weighted
   predictor second derivatives. Use the existing safeguarded mode/factor
   machinery and require undamped positive final curvature for Laplace.
3. Compute complete outer Laplace scores, including mode motion, determinant
   derivatives, cross-layer third-order contractions and moving covariance
   maps. Do not substitute first-order backpropagation or Gauss-Newton
   curvature and label it the existing exact-Hessian estimator.
4. Check independent dense tiny-network modes, Hessians and marginal
   likelihood/score references. Use numerical integration for very small
   random dimensions where practical to characterize Laplace error separately
   from derivative error.
5. Validate conditional, zero-deviation and joint-prior population prediction
   and simulation with checked integration error and explicit dimension/cost
   admission. Reject impractical integration instead of silently substituting
   zero effects or posterior-mode covariance.

**Exit gate:** Gaussian random-network recovery and complete derivative
references, declared covariance/RR coverage, deterministic mode behavior,
all point-output surfaces, serialization and bounded fit/integration memory.
No inference validity is inferred from successful optimization alone.

### Phase 11 — graph families, post-fit inference and release

**Prerequisite:** Phase 10. **Deliverable:** a documented graph support matrix
and only the statistically justified full user-facing surface.

1. Enable existing native scalar families, then categorical/mixed families,
   only after each passes complete objective/score/Hessian, nuisance-parameter,
   prediction and simulation checks. Preserve block likelihood curvature and
   existing trait/category conventions; do not claim automatic family coverage
   merely because forward evaluation works.
2. Reuse explicit post-fit uncertainty operations for identifiable parameter
   coordinates and prediction functionals. Propagate mode sensitivity,
   covariance/RR maps and preprocessing. Reject singular or symmetry-dependent
   weight covariance, unsupported intervals and uncalibrated comparisons;
   prediction uncertainty also requires its own estimability check.
3. Keep PEV, observed parameter curvature, delta-method covariance and intervals
   outside fitting. Distinguish conditional covariance from approximate PEV and
   report Laplace/delta-method limitations. Do not advertise ordinary Wald/LRT
   calibration for hidden-unit counts, boundary variances or unidentified weights.
4. Qualify scaling by depth, width, observations, groups and random parameter
   count; report derivative, mode, optimization and explicit post-fit costs
   separately. Full-batch likelihood, sparse random algebra and memory guards
   remain the baseline. Alternative optimizers, accelerators or approximations
   need separate approval and reference gates.
5. Run repeated-simulation recovery/coverage studies as standalone workloads,
   all classical/graph regressions, platform CI, package/manual/vignette checks,
   fresh-process serialization and API/capability documentation. Every automated
   case remains below 30 seconds.

**Exit gate:** all enabled graph/family combinations meet declared numerical,
statistical, output, resource and portability criteria. Unsupported surfaces
fail explicitly. Multi-layer support is complete only for this qualified
envelope, not arbitrary neural architectures.

#### Later research — not a committed phase

ReLU and other nonsmooth operations need a separate policy for knots and
curvature before random-parameter Laplace or uncertainty can be supported.
First-order subgradients alone do not supply that policy. Wider/deeper models,
penalization, stochastic training, accelerators and richer architectures are
separate research decisions; none bypasses the classical or graph gates.

## Main technical risks and mitigations

| Risk | Required mitigation |
| --- | --- |
| Non-convex fixed network likelihood and equivalent parameterizations | Document coefficient non-uniqueness; deterministic starts; expose convergence and gradients; add controlled multi-start only after cost and selection semantics are specified. |
| ReLU is non-differentiable at zero | Specify the zero convention; test kink behavior; report that local curvature/inference at a kink is unavailable or approximate. Include smooth activations for inference comparisons. |
| Nonlinear random-mode Hessian differs from current linear Hessian | Include predictor second derivatives in exact Hessian and Laplace determinant; do not substitute Gauss-Newton curvature without labeling the resulting approximation. |
| Outer score must include mode and determinant dependence | Validate full-objective scores against central differences. Use numerical scores as a reported initial route if analytic/automatic scores are not yet proved. |
| Random network dimensions can grow as units × groups | Start with small group-level parameter blocks, enforce memory checks, and reject unsupported dense layouts before allocation. |
| Standard Gaussian random-effects assumptions conflict with parameter transforms | Define prior on the exact parameter scale; document transforms and covariance interpretation. |
| Existing result methods assume formula-based linear designs | Retain a serializable predictor graph and add explicit nonlinear result/prediction branches; do not forge a linear `X` that implies false coefficient semantics. |
| Statistical model-selection counts are ambiguous for constrained/symmetric networks | Count free independent fitted coordinates consistently; defer AIC/BIC or tests if the count or likelihood is not defensible. |

## Recommended next implementation step

Complete Phase 7 classical-model consolidation, starting with the outstanding
external release gates and a full regression run against Phase 6. Agree the
support matrix and numerical/statistical thresholds before extending tests.
Resolve supported-scope defects before broadening model classes.

Only after that gate and explicit approval, proceed through the graph and
backpropagation foundation (Phase 8), fixed smooth multi-layer models (Phase 9),
small random network blocks (Phase 10), and qualified families/post-fit
inference/release (Phase 11). The optional classical backlog is prioritized
separately. This ordering preserves modularity, existing formulas, validated
control options and strictly post-fit PEV/uncertainty throughout.
