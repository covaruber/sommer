# Nonlinear mixed-model effects: phased implementation plan

## Progress checkpoint (2026-11-11)

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
| 7: Classical-model consolidation and release gate | Skipped by user direction; exit gate remains open | Portability, full-suite, numerical/statistical and release qualification are not claimed. |
| 8: Declarative graph and derivative foundation | Implemented as bounded standalone prototype by explicit user direction | Versioned serializable graph evaluator, fixed/random parameter maps, native reverse gradients, Hessian-vector products and third-order directional contractions; no network fitting. |
| 9: Fixed multi-layer predictors | Implemented and locally validated; Phase 7 remains skipped and unqualified | Fixed smooth graphs, exact Gaussian ML, conventional additive Gaussian random effects, prediction and simulation. |
| 10: Random multi-layer parameters | Implemented and locally validated for the stated Gaussian scope | Independent and joint graph-coordinate covariance blocks, exact conditional curvature, numerical outer Laplace scores, prediction, simulation and reduced-rank maps. |
| 11: Unified specifications, graph families, post-fit inference and release | Implemented and locally qualified within the stated envelope | Native/reference scalar, categorical and mixed families, outputs, checked local inference and resource admission; external release gates remain open. |
| 12: Release baseline and classical stability closure | Planned; not implemented | Close Phase 7 and outstanding platform/manual gates before expanding the advertised release envelope. |
| 13: ReLU and other nonsmooth activations | Planned; not implemented | Explicit knot semantics, smooth comparators, and separate fitting/integration/inference gates. |
| 14: Wider networks and scalable graph execution | Planned; not implemented | Extend measured capacity beyond current bounds through derivative/solve improvements, not by raising limits alone. |
| 15: Penalized multi-layer mixed models | Planned; not implemented | Explicit fixed-weight penalties, tuning, optimizer and uncertainty semantics distinct from ML and random-effect priors. |
| 16: Stochastic training | Planned; not implemented | Statistically justified sampling units, reproducible stochastic optimization and full-target certification. |
| 17: Graph accelerators | Planned; not implemented | Optional operation-local device kernels with CPU parity, precision and memory admission. |
| 18: Extended-model qualification and publication | Planned; not implemented | Cross-feature statistical, numerical, resource and platform evidence before advertising combined support. |

Roadmap entries are not automatic authorization for later implementations.
On 2026-11-09, the user directed that Phase 7 be skipped and Phase 8 proceed.
The user later authorized full Phase 9 and Phase 10 implementations and
acknowledged that Phase 7 remains skipped and unqualified. Phase 11 was
subsequently authorized and implemented; it does not close Phase 7's gates.
On 2026-11-11 the user requested the roadmap for Phases 12--18 and explanatory
multi-layer vignettes, not implementation of those new engine features.

Phase 11 was revised and implemented against the unified public interface.
Earlier phase examples and validation counts describe their
historical engine contracts; they are not current public call syntax or
evidence that the newer specification route has the same qualifications.
Current calls use formula expressions or `graph_term()` with `parameters`
and `predictors`, `method=mmes_method()`, top-level natural-scale `start`,
and `control=mmes_control()`. The former public `nonlinear=` argument and
`control$task` module have been replaced, not retained as aliases.

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
The Phase 7 prerequisite was explicitly waived for this bounded prototype;
the open classical stability risks remain and do not transfer to Phase 9.
**Deliverable:** a bounded native graph/derivative prototype, not a new
production estimator.

Implementation and local qualification:

- Added public `nonlinear_graph()` and `nonlinear_graph_eval()` functions.
  Graphs are versioned, serializable DAGs with declared input transforms,
  affine nodes, identity/tanh activations, shared parameter coordinates, and
  fixed/group-random maps. Evaluation does not fit parameters.
- Native evaluation uses bounded workspace admission, reverse-mode gradients,
  Jacobian-vector products, Hessian-vector products, third-order directional
  gradients and directional derivatives. Directional operations avoid dense
  third-order parameter tensors.
- Identity-graph/linear and tanh-composed logistic parity passed. Analytic
  finite-difference checks cover multilayer paths, shared/transformed
  coordinates, group-specific directions and third-order contractions.
  Serialization passed a fresh-process check.
- The dedicated graph tests passed 32 assertions. All nonlinear regression
  tests passed 2,626 assertions with no failures, warnings or skips on the
  local configuration. This does not qualify other platforms.

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
bounded memory; unchanged classical regression suite. **Prototype exit gate
passes locally.** Phase 7 remains skipped and unqualified. No network fitting
is advertised until Phase 9 passes.

### Phase 9 — smooth fixed multi-layer predictors

**Status:** implemented and locally validated within the bounded fixed-network
scope below. Phase 7 remains skipped and unqualified; Phase 9 does not close
its release gates.

**Prerequisite:** Phase 8 graph evaluator. **Deliverable:** opt-in exact Gaussian
ML fitting for a small smooth fixed network plus conventional additive Gaussian
random effects. Existing nonlinear-curve routes and optimizer are preserved.

#### Implemented public and statistical contract

- Accept a `sommer_nonlinear_graph` through existing nonlinear-model fitting
  dispatch. Keep `nonlinear_curve()` behavior and dispatch unchanged.
- In this phase, require every graph parameter coordinate to be fixed.
  Reject random graph coordinates explicitly; defer random network parameters
  to Phase 10.
- Start with `tanh` hidden activations and identity output, one hidden layer,
  and a bounded two-hidden-layer acceptance graph. Reuse current declared
  graph operations only; do not add arbitrary callbacks or operations.
- Combine graph output with the ordinary fixed design, offset, and additive
  Gaussian random effects. Integrate additive random effects exactly and
  optimize fixed network, conventional fixed-effect, covariance, and residual
  variance parameters under Gaussian ML.
- Require a complete named initial vector for graph parameters in this first
  fitting version. Keep initialization deterministic; defer random
  initialization, multistart, and new control fields unless implementation
  evidence shows they are required.
- Use graph-declared input transforms unchanged at fit and prediction time.
  Do not learn preprocessing parameters implicitly.
- Do not compute fit-time PEV or parameter covariance. Do not claim global
  optimum, unique network weights, or inference validity from convergence.

#### Implemented components

1. **Graph fitting preparation.** Extended the current nonlinear fit
   dispatcher to distinguish curve models from graph objects. Reuse existing
   formula preparation, missing-data handling, fixed/random design creation,
   covariance descriptors, and workspace admission. Validate graph input
   columns, fixed-only parameter maps, finite named starts, response shape,
   and graph/fixed-design identifiability before optimization.
2. **Exact Gaussian ML.** Connected graph output and analytic Jacobians to the
   existing native Gaussian objective. Reused existing Gaussian
   covariance integration and covariance/residual scores; do not duplicate
   likelihood algebra or substitute squared-error training, PQL, or Laplace
   integration for the exact Gaussian route. Kept `.nonlinear_optimize()` and
   `nlminb`.
3. **Independent objective and score proof.** Added tiny dense reference
   calculations for the marginalized Gaussian objective, including a
   conventional random intercept and nontrivial covariance. Verified native
   scores against independent finite differences for fixed, graph, covariance
   and residual-variance coordinates, including cross-layer paths and transforms.
4. **Special-case and fit acceptance.** Verified identity-graph semantics and
   fit parity for a graph-composed classical logistic curve. Added one- and
   two-hidden-layer simulated fixtures; assess recovery by predictions or
   identifiable functionals, not raw weights.
5. **Identifiability and convergence.** Check the augmented fixed
   design/graph Jacobian at initialization and final fit. Detect local rank
   deficiency and confounding with ordinary fixed terms. Document hidden-unit
   permutation and scaling symmetries; do not canonicalize equivalent
   networks or reject them using arbitrary weight ordering. Report optimizer
   convergence, objective, gradient, starting values, and existing rejected
   evaluation diagnostics.
6. **Fitted-object behavior.** Persist the serializable graph,
   fitted parameter coordinates, input-transform metadata, formulas/design
   metadata, and standard fit diagnostics; never persist native pointers.
   Implemented coefficients, summary/print, fitted values, residuals,
   likelihood, new-data prediction, simulation, and profiling through existing
   nonlinear methods. Graph-parameter uncertainty requests fail explicitly.
7. **Resource and regression gates.** Bound graph, observation,
   parameter, tape, derivative, and prediction workspace against
   `execution$denseMemoryMB`; excess fails explicitly. Tested deterministic
   repeated fits and unchanged caller RNG state. The complete nonlinear
   regression suite passes without changing classical curve behavior.

**Exit gate:** linear and classical-curve parity; independently verified
marginal Gaussian objective and all score blocks; one- and two-hidden-layer
fit acceptance; identifiable prediction/function recovery; deterministic
repeat and RNG preservation; serialization and new-data prediction parity;
complete fitted-object surface; memory rejection; and unchanged nonlinear
regressions. Local Phase 9 success would not close Phase 7 release gates.

**Local validation:** graph evaluator and fitting acceptance passed 79
assertions. The complete nonlinear regression suite passed 2,673 assertions
with no failures, warnings, or skips on the local configuration. Cross-platform
qualification remains open with Phase 7.

**Explicit non-goals:** random network parameters, Phase 10 curvature and
Laplace derivatives, non-Gaussian families, ReLU, dropout, batch normalization,
stochastic/minibatch training, GPU execution, automatic multistart, global
optimality, and fit-time parameter uncertainty.

### Phase 10 — small random multi-layer parameter blocks

**Status:** Implemented and locally validated for the supported envelope.
Graph fits may opt in to Gaussian identity-link ML with explicitly declared
group-level random deviations on smooth parameters. Independent coordinates
use named `random` formulas; cross-parameter blocks use `joint` with a
`.parameter` covariance axis. Existing `vsm()` covariance, relationship and
reduced-rank machinery is reused. Priors are on graph working coordinates,
including log-transformed parameters.

1. The conditional mode uses the exact Hessian, including score-weighted
   predictor second derivatives, safeguarded optimization and a positive
   undamped final curvature requirement.
2. The outer Laplace score is central finite differences of the exact fitted
   objective, not an analytic score. This is recorded in fit diagnostics and
   remains a limitation to address before broader performance/release claims.
3. Conditional and population prediction, conditional and population
   simulation, serialization, independent/joint covariance blocks and reduced
   rank are covered. Population integration has explicit memory/order gates
   and is checked against an independent high-order one-coordinate reference.
4. Dense numerical references check the conditional mode, Hessian and
   marginal likelihood. The graph regression tests also check recovery on a
   seeded small grouped example and exact-score finite-difference agreement.
5. Graph non-Gaussian families and post-fit parameter uncertainty remain
   explicitly unsupported; unsupported requests fail rather than falling back
   to zero deviations or linear-model machinery.

The graph regression gate passes, and all nonlinear regression files pass
locally (2,705 assertions). Phase 7 portability, full-package-suite and
external release qualification remain open; successful Phase 10 fitting is
not evidence of general inference validity.

### Phase 11 — unified specifications, graph families, post-fit inference and release

**Status:** Implemented within the supported envelope below. Local statistical
and regression evidence is recorded here; broader release/platform gates
remain explicitly unqualified. The deliverable is one documented statistical
contract with explicit support by fitting route, family, covariance and output.

#### Implemented scope

| Workstream | Implementation |
|---|---|
| 11A: public contract | Preserved ML/REML/PQL routing and numerical-control ownership; transported public family nuisance declarations through both routes; retained ordinary Gaussian `fixedRank="drop"`. |
| 11B: representation | Serializable parameter/graph layouts, natural versus coefficient starts, declared transforms, retained-row/category mappings and route parity; reference RR remains explicitly rejected. |
| 11C: scalar families | Shared Gaussian, three binomial links with grouped trials, Poisson, Gamma, NB2 and ZIP kernels; exact nonlinear conditional curvature and reported central outer scores. |
| 11D: categorical/mixed families | Ordered ordinal thresholds, reference-category nominal expansion, full cross-logit curvature and mixed-family output layouts. |
| 11E: outputs | Conditional/zero/population prediction and simulation, joint predictor-prior quadrature without a two-record cutoff, coefficient/residual adapters and seed restoration. |
| 11F: uncertainty | Reconstructed and checked observed information, coefficient/natural/group covariance, mode sensitivity, approximate empirical Bayes PEV, intervals and guarded regular comparisons; no pseudoinverse repair for unidentified weights. |
| 11G: resources/documentation | Route-aware performance profiling, admission before dense/quad allocations, persistent Unix/Windows header dependencies, help, articles, capability metadata and generated decision guide. |
| 11H: qualification | Public-route, independent curvature/root/quadrature, serialization and recovery/coverage workloads; local evidence and outstanding release gates are distinguished below. |

Reference conditional derivatives now compose symbolic expression derivatives,
native graph derivatives and full native family curvature. Predictor second
derivatives and parameter-link second derivatives are retained. The outer
score remains numerical; this does not make a Laplace integral exact.
Legacy noncanonical scalar links retain a separately identified numerical
inverse-link second derivative.

Post-fit differentiation reconstructs conditional modes at tighter precision
than ordinary fitting, retaining the saved fitting controls and statistical
model. This fixes observed-information instability from differencing
insufficiently precise modes without relaxing derivative-agreement checks.
No uncertainty is computed during fitting.

Local evidence on Linux/R 4.2.2/FlexiBLAS includes:

- 182 working-tree nonlinear/public/specification cases: 3,131 expectations,
  no failures, errors, warnings or skips.
- Full working-tree suite: 836 cases, 17,198 expectations, no failures or
  errors. One existing native-Laplace false-convergence warning remains.
  The two skips are GPU qualification and the installed-namespace-only gate;
  that namespace gate subsequently passed separately against the freshly
  checked installation. Numerical validation used process-local single-thread
  BLAS/OpenMP settings rather than changing shared-system configuration.
- The RR acceptance fixture retains four groups, independent/joint mapping
  and restored-versus-latent checks; its measured case is 16.952 seconds.
  The complete graph-random file remains below 30 seconds per case.
- Independent scalar-root binomial Laplace information references at seeds
  75110 and 75123, nonlinear graph conditional derivatives checked against
  finite differences, and a one-dimensional quadrature bound on Laplace error.
- Fresh-process native/reference prediction, covariance, confidence intervals
  and seeded simulation reconstruction.
- Source-package installation, help, S3 consistency, compiled-code checks and
  examples passed. The full vignette build, execution and output rebuild
  passed separately. After qualifying `stats::deriv3`/`stats::vcov`, the final
  source check (`--no-manual --no-tests --ignore-vignettes`) has no errors or
  warnings and one installed-size NOTE; its R-code usage check is clean.
  Tests and vignettes are qualified by the separate runs above, not by those
  skip flags.
- Generated guide checks passed: 115 options, 18 nodes, 15 rules, 147 fit-time
  checks, 18 structure paths and 12 tips. Both affected articles and all changed
  help sources were rebuilt/checked.
- `inst/benchmarks/nonlinear-phase11-qualification.R`: 30 replications for
  Gaussian, binomial and Poisson on each route, 180 successful fits/inference
  requests with no failures. The estimand is the zero-deviation link predictor
  at `x=0`, with prior random-intercept variance `.09` known.

| Family | Native/reference bias | Native/reference 95% coverage |
|---|---:|---:|
| Gaussian | -0.01819079 | 0.9666667 |
| Binomial | 0.11059059 / 0.11059058 | 0.9333333 |
| Poisson | -0.03486263 | 0.9333333 |

These are finite Monte Carlo checks for regular identified estimands with
known covariance, not general coverage claims for arbitrary network weights,
estimated boundary variances or all categorical interval functionals.
Object sizes in the workload are not peak resident memory.
Eight unmodified legacy full-suite cases still exceed 30 seconds; the new and
edited Phase 11 cases meet the measured per-case limit. These local runs do not
close older global runtime or platform gates.

The PDF manual is not qualified locally: `pdflatex` is absent, and a private
TinyTeX installation attempt failed at the HTTPS download. Windows, macOS,
R-devel, sanitizers and other external/platform release gates remain unexecuted.
Phase 7 remains skipped and unqualified.

#### Pre-implementation baseline and changes to the earlier plan

The public entry point is `R/FUN_model_interface.R`. It selects native
components or the bounded reference engine in `R/FUN_model_specification.R`.
These are not interchangeable implementations of every model:

| Route | Behavior before Phase 11 |
|---|---|
| Ordinary Gaussian models, including fixed-only | Baseline engine; fixed-rank handling and covariance output are retained. The earlier fixed-only Gaussian dispatch regression has a dedicated public-interface regression test. |
| One additive native component with no parameter overrides | `graph_term("net")` and `predictors=list(net=net)` can select the existing native curve/graph engine. Native graph internal random maps, independent formulas and `joint` blocks remain available on this route. |
| Fixed-only native component without random effects | Redirected to the reference specification route; a native graph declaration alone does not guarantee a native fitting engine. |
| Expressions, parameter overrides, or other graph compositions | Reference specification route, with namespaced coordinates such as `"net::k"=parameter(fixed=~Env, random=~group, link="log")`. Graphs must have no internal random maps or `random`/`joint` declarations. |
| Reference engine limits | Numerical outer/mode scores and conditional Hessians, dense bounded workspaces, scalar Gaussian/Bernoulli/Poisson/Gamma likelihoods; graph fitting remains Gaussian identity only. Nonempty `familyParameters` and many non-default controls are rejected. |
| Result methods | Native nonlinear and specification results have different classes and adapters. Reference population prediction currently admits at most two random predictor records; simulation is explicitly unavailable. Reference `vcov()` computes post-fit numerical curvature but is not a qualification of general graph inference. |

Consequently, enabling a native graph family alone would have left common
`parameter()` calls unsupported. Conversely, removing the reference engine's
Gaussian graph guard would not establish native derivative, covariance,
resource or inference parity. Both routes require explicit gates.

The statistical target is
`g(E[Y_i | b]) = X_i beta + h_i({T_j(A_ij alpha_j + Z_ij b_j)}) + Z_0i b_0`,
where `h_i` can include supported expressions and graph components, `T_j`
is the declared identity or log parameter link, and the latent effects have
their declared Gaussian covariance. Parameter links and the response link
are distinct. Separate parameter random records are not implicitly a joint
cross-parameter covariance block merely because they share a grouping factor.

#### 11A. Freeze public semantics, dispatch and control ownership

1. Build a support matrix indexed by ordinary/native/reference route,
   expression versus graph, fixed-only/additive/parameter-random structure,
   family/link, covariance/RR map, and requested output. Include execution
   limits and whether scores/conditional curvature are analytic or numerical.
2. Keep statistical choices in `mmes_method()`, nuisance declarations in
   `familyParameters`, starts in `start`, and numerical choices in the public
   optimizer/mode/likelihood/population/inference control modules. Do not add
   another public fitter or resurrect `nonlinear=` or `control$task`.
3. Preserve Gaussian REML and mixed non-Gaussian PQL defaults. Nonlinear
   likelihood fits require explicit ML; non-Gaussian random models cannot
   request exact integration. Fixed-only likelihoods have no random integral;
   Gaussian identity models linear in their Gaussian effects retain exact
   integration even when their fixed mean is nonlinear.
4. Report the selected route, actual integration and derivative methods,
   dimensions and admission reason. Equivalent supported declarations must
   not silently change the statistical model or likelihood normalization.
5. Test controls at the selected route, not merely at construction. Preserve
   `representation$fixedRank="drop"` for ordinary Gaussian models. For
   nonlinear models, distinguish redundant linear columns from genuine
   nonlinear non-identifiability; do not drop arbitrary network coordinates
   or remove blanket control guards without implementing their semantics.
6. Preserve both existing random-declaration conventions without permitting
   double specification. Decide any new reference-route joint-parameter API
   before implementation; do not invent shared covariance from matching names.

**Gate:** public dispatch/control tests cover each row of the support matrix,
explicit unsupported requests, and unchanged ordinary Gaussian/PQL behavior.

#### 11B. Qualify preparation, parameter maps and Gaussian route parity

1. Record a serializable model layout containing graph namespaces, repeated
   coordinates, parameter fixed regressions, random design records, offsets,
   covariance coordinates, transforms and free family coordinates. Reuse
   existing layouts where possible rather than introducing a new general IR.
2. Preserve formula environments, contrasts, levels, relationship precision
   ordering, one retained-observation mask, preprocessing and prediction
   templates. Apply parameter links exactly once: for log parameters,
   `theta_i=exp(A_i alpha+Z_i b)`, with Gaussian priors on the working scale.
3. Validate natural parameter starts separately from named coefficient starts,
   including non-intercept regressions and qualified graph names. A regression
   coefficient start is not a vector of natural graph parameter values.
4. Compare native and reference formulations only where they describe the
   same model: fixed-only, exact Gaussian additive effects, and small
   independent random-coordinate Laplace models. Check means, covariance,
   likelihood constants, modes and target predictions, not just convergence.
5. Qualify relationship and reduced-rank maps per route. The dense reference
   covariance/precision path is not automatically suitable for a singular RR
   covariance; reuse proper latent coordinates or reject that combination.
   Do not add a ridge to manufacture an invertible prior.

**Gate:** Gaussian public-interface parity and transform/namespace tests pass,
with route-specific unsupported covariance cases stated explicitly.

#### 11C. Enable scalar graph families through both public routes

1. First qualify fixed-only binomial-logit and Poisson-log graphs on the
   reference route. Then qualify additive Gaussian random effects on native
   and reference routes, followed by random graph coordinates and parameter
   regressions. Use the same small scientific predictor wherever routes overlap.
2. Reuse native family response validation, likelihood constants, score,
   curvature and nuisance contracts rather than maintaining a second family
   registry. Adapt reference preparation and result assembly as well as native
   graph dispatch; do not only remove the Gaussian eligibility guards.
3. Add binomial probit/cloglog, Gamma-log, NB2-log and ZIP-log in individually
   gated increments. Grouped binomial trials are separate from the reference
   engine's current Bernoulli-only contract. Fixed/estimated dispersion,
   NB2 shape and ZIP structural-zero predictors need explicit layout tests.
4. Native conditional curvature must include the full observation curvature
   and score-weighted predictor second derivatives. Check parameter-regression
   and graph-transform chain rules, mixed derivatives, mode motion and moving
   covariance/RR maps. Require positive undamped curvature at the accepted mode.
5. Retain central finite differences as a clearly reported outer-score method
   where necessary. Reference conditional Hessians compose symbolic/native
   chain rules and are checked against numerical mode-score derivatives.
   The legacy scalar-link numerical exception is identified separately. Analytic/native
   scores are a separate deliverable requiring third-order and nuisance
   references; never substitute Gauss-Newton for the exact conditional target.
6. Validate likelihood, score and Hessian against independent tiny references,
   including nonconcave ZIP starts. Use low-dimensional numerical integration
   to characterize Laplace error separately from derivative error.

**Gate per family:** fixed-only, additive and parameter-random cases pass
their declared numerical, output and resource checks before being advertised.
Native success does not enable reference support, or vice versa.

#### 11D. Qualify categorical and mixed-family representations

1. Add ordinal logit/probit only with explicit threshold ordering, free-location
   identifiability checks, response coding and nuisance derivatives.
2. Design nominal category-specific predictor rows before enabling multinomial
   graphs. A scalar graph can supply a scalar value per expanded logit row;
   category-specific fixed/random mappings, shared parameters, reference
   category and cross-logit curvature still need an explicit design.
3. Generalize the reference scalar-response contract only after specifying
   categorical/mixed response layout and outputs. Reuse existing native
   `familym` trait/category indexing; do not infer support from forward graphs.
4. Preserve full multinomial block curvature and cross-category random
   covariance in mode fitting, scores, population integration and simulation.

**Gate:** probability sums, category/reference ordering, thresholds, missing
cells and mixed-family parameter mappings pass independent references.
Categorical support remains a Phase 11 deliverable, not silently removed from
scope; unsupported combinations remain disabled until these gates pass.

#### 11E. Complete route-aware prediction, simulation and result adapters

1. Preserve public `target="conditional"`, `"zero"` and `"population"` and
   link/response output semantics for both result classes. Population link
   output means the prior expectation of the nonlinear link-scale predictor,
   not necessarily its zero-deviation value.
2. Integrate the actual joint prior of the low-dimensional predictor coordinates,
   including graph parameters, additive effects, correlated parameter blocks
   and categorical logits. Reuse covariance/design maps and numerical error
   checks. Replace the reference two-record limit only with validated
   dimension/cost admission; do not integrate a high-dimensional group vector.
3. Implement reference-route conditional and population simulation, currently
   rejected. Draw the fitted prior covariance, apply parameter links and graph
   transforms once, then simulate the family response. Check seed restoration,
   category/trial conventions and empirical moments.
4. Make coefficient labels, fitted values, residuals, summaries, likelihood
   degrees of freedom and profiling consistent without forcing incompatible
   result classes into curve-only fields. Test fresh-process serialization;
   never persist native pointers or factors.

**Gate:** each enabled route/family has complete point outputs and checked
population integration; unsupported targets fail rather than returning modes
or zero effects as a success-shaped substitute.

#### 11F. Add explicitly post-fit inference for identified estimands

1. Reconstruct either route's saved likelihood and verify its value before
   inference. Extend shared post-fit evaluation adapters to specification
   results; do not simply remove native graph rejection checks.
2. Qualify the existing reference `vcov()` path against the same convergence,
   score, boundary, competing-mode, derivative-stability and positive-curvature
   checks as other inference. Do not invert a raw numerical Hessian and infer
   statistical validity solely from successful factorization.
3. First support regular, locally identified fixed/parameter-regression
   coordinates and prediction functionals. Include parameter and response-link
   transforms, estimated nuisance parameters, preprocessing and covariance maps.
   Reject singular or symmetry-dependent raw network-weight covariance.
4. Prediction inference requires an estimability check for the requested
   functional. Do not use a generic pseudoinverse to disguise unidentified
   nuisance directions; any estimable-subspace method needs its own reference
   and parameterization-invariance tests before admission.
5. Add conditional random-mode covariance and approximate empirical Bayes PEV
   with mode sensitivity, cross-covariance, RR restoration and group-parameter
   regression derivatives. Distinguish inverse conditional curvature from
   parameter-adjusted PEV and from exact frequentist prediction-error variance.
6. Propagate population integration and estimated-parameter uncertainty for
   each interval target. Keep all SE/PEV/interval computations outside fitting.
   Restrict comparisons to justified regular cases; no ordinary Wald/LRT
   calibration for graph width, boundary variances or unidentified weights.

**Gate:** dense curvature/mode-sensitivity/delta references and standalone
coverage studies pass for each advertised estimand. Unsupported uncertainty
fails explicitly in both result classes.

#### 11G. Qualify resources and publish the route-specific contract

1. Measure preparation, graph derivatives, modes, outer scores and post-fit
   work separately by observations, groups, depth/width, parameter-regression
   coefficient count, latent rank and family. Include the reference engine's
   current dense limits (30 declared parameters, 60 mean coefficients and
   512 modes for non-exact models) rather than claiming native sparse scaling.
2. Keep admission before dense allocations and quadrature expansion. Native
   sparse factor/solve machinery remains the scalable baseline; retain a bounded
   reference path for verification. Do not raise limits to imply scalability.
   Any native backend extension must respect signed/block curvature and
   operation-local solve adapters, not assume a stored random-effect Hessian.
3. Update `R/FUN_mme_capabilities.R`, the numerical controls reference/vignette,
   `vignettes/sommer.models.Rmd`, graph/parameter/method help, and the generated
   options guide together. Show accepted values separately from availability
   by route, including fixed-rank and likelihood-backend controls.
4. Keep alternative optimizers, stochastic training, GPU execution, new graph
   builders and nonsmooth activations outside this implementation scope.

**Gate:** measured resource limits and runtime checks agree with all published
support tables; no ignored non-default controls or undocumented fallbacks.

#### 11H. Regression, statistical and release qualification

1. Extend existing `test-model-interface.R`, `test-model-specification.R` and
   graph/family/inference tests with a small public-route acceptance matrix.
   Keep private numerical-engine fixtures: they test different boundaries.
2. Use independent numerical oracles and simulated truth without adding
   `lme4`/`nlme` dependencies. Require mode/outer-score certificates, objective
   parity and declared error tolerances rather than optimizer return codes alone.
3. Keep every automated case below 30 seconds. Put repeated-simulation
   recovery/coverage and scaling in standalone `inst/benchmarks/` workloads,
   recording seeds, failures, approximation error and interval targets.
4. Run working-tree focused tests, all nonlinear/public-interface regressions,
   then the full package suite, source-package/manual/vignette checks,
   fresh-process serialization and supported-platform CI. Re-establish current
   evidence; historical Phase 10 assertion counts do not validate new dispatch.

**Dependency order:** 11A -> 11B -> 11C -> 11D; 11E follows each enabled family,
11F follows its validated likelihood and point outputs, and 11G runs alongside
each increment. 11H release qualification depends on all enabled combinations.
Do not defer basic output or documentation gates until the final family.

**Exit gate:** all advertised route/family/covariance/output combinations meet
their declared numerical, statistical and resource criteria. Record unexecuted
platform gates explicitly. Phase 7 remains skipped and unqualified: local
Phase 11 success does not retroactively close its release requirements.
Multi-layer support is complete only within the qualified envelope, not for
arbitrary neural architectures or arbitrary formula callbacks.

### Extension terminology and scope

Use **multi-layer nonlinear mixed model** as the umbrella term when a
multi-layer predictor has random effects. For non-Gaussian families, use
**generalized nonlinear mixed model with a multi-layer predictor**. Without
random effects these are nonlinear/generalized nonlinear regressions, not
mixed models. Layers themselves are deterministic; a Gaussian random-weight
model requires an explicit covariance declaration.

Smooth multi-layer graphs already exist. The phases below extend their
activation, capacity, objective or execution envelope rather than introduce
multi-layer models for the first time. They preserve the formula,
`parameters`/`predictors`, `method`, `start` and `control` separation.
New syntax must be designed and approved at its phase boundary; the roadmap
does not promise currently nonexistent controls or convenience builders.

### Phase 12 — release baseline and classical stability closure

**Status:** planned, not implemented. **Prerequisite:** the Phase 11 local
baseline. **Deliverable:** a reproducible release-qualified reference
envelope against which the new extensions are compared.

1. Reopen and complete the skipped Phase 7 checklist against current sources,
   retaining historical validation rather than relabeling skipped work complete.
2. Run supported Windows/macOS/Linux and R-release/R-devel checks, sanitizer
   builds where supported, installed-namespace checks, fresh-process outputs,
   package/vignette builds and the PDF manual on a TeX-equipped host.
3. Inventory the existing warning/skips and legacy tests exceeding 30 seconds.
   Resolve defects or document justified qualifications; preserve their
   numerical coverage when moving scaling work to standalone benchmarks.
4. Publish one support/evidence matrix crossing route, graph depth/width,
   family/link, random-coordinate structure, covariance/relationship map,
   output, derivative method and platform. Untested intersections stay unqualified.
5. Capture deterministic CPU likelihood/score/mode/inference references,
   timings, peak memory and toolchain metadata for use by later phases.

**Gate:** reproducible release evidence and an explicit disposition for every
Phase 7/platform/manual item. Local successes alone do not close this gate.

### Phase 13 — ReLU and other nonsmooth activations

**Status:** planned, not implemented. **Prerequisite:** Phase 12.
**Deliverable:** individually admitted activation/model combinations with
honest differentiability, integration and uncertainty semantics.

1. Define versioned node contracts for ReLU, then selected additional
   activations only after separate validation. Specify forward values,
   derivative/subgradient conventions at zero, parameter domains and
   serialization. A smooth comparator such as softplus must be named and
   parameterized as a different model, never substituted silently for ReLU.
2. Start with fixed graph coordinates and no random effects or additive
   Gaussian effects, where observation-space activation regions can be
   examined independently. Implement an optimizer suitable for piecewise
   smooth objectives, with one-sided checks and an explicit generalized
   stationarity certificate at knots; ordinary smooth-gradient stopping is
   insufficient there.
3. Treat random coordinates crossing activation boundaries as a separate
   integration problem. Study both interior-region and knot modes against
   one-dimensional adaptive quadrature and partitioned tiny references,
   including multiple activation regions. A piecewise Hessian away from a
   knot does not by itself validate a Laplace integral across regions.
4. Decide and document the admitted treatment of knot modes: justified
   region-aware integration or an explicit unsupported request. Ordinary
   smooth Laplace, damped positive curvature and first-order subgradients
   are not interchangeable solutions. Keep response-family signed/block
   curvature and proper latent covariance coordinates intact.
5. Gate uncertainty independently from point fitting. Smooth local
   information is available only under justified stable-active-region
   assumptions. At knots or active-set changes, reject ordinary Wald
   inference unless a separately validated nonregular procedure is provided.
6. Propagate activation semantics through prediction, population quadrature,
   simulation, profiling and save/reload; add knot-adjacent and boundary
   fixtures that fail if an activation or its derivative is silently smoothed.

**Gate:** each advertised activation has independent value/derivative,
optimization, integration and output references. Unsupported random-knot or
uncertainty cases remain explicit errors, not a phase-wide claim of arbitrary
nonsmooth mixed-model support.

### Phase 14 — wider networks and scalable graph execution

**Status:** planned, not implemented. **Prerequisite:** Phase 12; Phase 13
only for nonsmooth variants. **Deliverable:** a measured larger smooth-graph
envelope with unchanged statistical targets.

1. Benchmark observation count, layer depth, layer width, tied/free
   coordinates, random-coordinate groups, relationship connectivity and
   categorical expansion separately. Record peak resident memory, derivative
   workspace, factorizations/solves, score/mode accuracy and whole-fit time.
2. Replace avoidable full derivative tensors with blocked/Jacobian-vector,
   Hessian-vector and directional contractions where the algorithm permits.
   Specify operator requirements before choosing dense versus sparse
   storage; retain predictor second derivatives and multinomial blocks.
3. Extend native operation-local factor/solve and latent RR mappings rather
   than assuming a stored random-effect Hessian or inverse. Preserve the
   bounded dense reference route as an oracle; increasing its numeric caps
   alone is not a scaling implementation.
4. Validate whole-objective outer derivatives, including mode and determinant
   dependence. If analytic/automatic derivatives replace central differences,
   qualify higher-order contractions, family nuisance terms and moving
   covariance maps against independent tiny references before switching defaults.
5. Add explicit initialization/scaling diagnostics for saturation and weak
   identification. If multi-start is added, expose starts, budget, selection
   criterion and every failed attempt; no opaque random restarts. A higher
   parameter count is not evidence that raw weights are identifiable.
6. Re-admit dense post-fit covariance, predictor-prior quadrature and
   simulation by their own dimensions/costs. Wider fitting support must not
   implicitly promise full weight covariance or high-dimensional population
   integration. Publish measured thresholds and fail before oversized allocation.

**Gate:** small-model likelihood/mode/output parity plus reproducible scaling
and memory evidence beyond the former envelope. Every raised limit is backed
by successful measurements and failure/admission tests.

### Phase 15 — penalized multi-layer mixed models

**Status:** planned, not implemented. **Prerequisite:** Phase 12 and a
declared resource envelope from Phase 14; Phase 13's nonregular policy for
activation/penalty combinations needing it. **Deliverable:** an explicit
penalized objective, not penalties disguised as Gaussian ML.

1. Specify the target
   `Q(lambda; rho) = -log L(lambda) + P_rho(alpha)` with the marginal likelihood
   and its constants intact. Define whether coefficients are penalized on
   natural or working scales, how tied coordinates are counted, which
   weights/biases are exempt and whether loss is summed or averaged.
2. Start with fixed-coordinate L2 penalties and fixed tuning values;
   include zero-penalty parity with current ML. Add L1/elastic-net or grouped
   penalties only with proximal/active-set methods and appropriate KKT
   certificates. Generic smooth BFGS at an L1 kink is not sufficient.
3. Keep fixed-weight regularization distinct from the existing Gaussian
   random-effect density and its normalization. Do not double-penalize
   random modes, conflate a tuning parameter with an estimated variance,
   or penalize covariance/nuisance coordinates without a separate contract.
4. Design grouped/nested cross-validation or a justified alternative for
   tuning. Fit preprocessing and select starts/architecture inside training
   folds. Hold out groups or connected units appropriate to the target;
   row-wise splits do not validate new-group prediction under dependence.
5. Report unpenalized likelihood, penalty, optimized criterion, tuning rule,
   active set and convergence separately. Do not label the criterion REML,
   ordinary `logLik`, AIC or a calibrated LRT without a justified definition
   of likelihood degrees of freedom and model-selection uncertainty.
6. Define uncertainty by estimand and interpretation. Penalized inverse
   curvature, frequentist sampling covariance and a MAP/posterior
   approximation are different objects. Qualify any sandwich/bootstrap or
   selection-aware procedure independently; otherwise reject unsupported
   interval requests rather than reusing unpenalized Wald calculations.

**Gate:** independent penalized tiny-model oracles, KKT/gradient checks,
zero-penalty regression parity, fold-leakage tests, serialized penalty/tuning
metadata and recovery/prediction studies reporting regularization bias.

### Phase 16 — stochastic training

**Status:** planned, not implemented. **Prerequisite:** Phase 14, plus
Phase 15 for penalized targets and Phase 13 for nonsmooth targets.
**Deliverable:** opt-in stochastic optimization with a stated statistical
target and reproducible approximation/certification policy.

1. Decide the sampling unit before choosing an optimizer. Row-wise unbiased
   gradients can be justified for factorized fixed-only likelihoods, but
   not generally for integrated mixed-model likelihoods. Start with
   independent complete groups/components when their marginal contributions
   factorize. Known relationships or shared latent effects can connect
   groups and invalidate that decomposition.
2. Derive and check gradient scaling for each admitted sampler/target.
   Re-solving a mode on an arbitrary subset and taking its log determinant
   is not an unbiased full Laplace score. Any approximate variational,
   subsampled or other objective must be named separately and never returned
   as exact ML/REML or the existing Laplace likelihood.
3. Introduce explicit optimizer, batch/sampler, learning-rate schedule,
   RNG stream and work budgets only after the target is defined. Preserve
   full-batch deterministic fitting as the default. Distinguish stochastic
   optimization from stochastic trace/linear-solve approximations already
   used by other engines.
4. Checkpoint coordinates, optimizer moments, sampler/RNG state, preprocessing,
   graph/family layouts and data-order compatibility, not native pointers.
   Test interrupted/resumed equivalence and avoid changes to the caller's RNG.
5. Monitor a reproducible full-data or independently justified diagnostic
   objective, score/KKT residual, mode feasibility and approximation error.
   A noisy training-loss plateau is not a convergence certificate. Where
   feasible, finish with full-objective deterministic refinement; if full
   certification cannot be computed, report training/approximation status,
   not certified marginal convergence.
6. Keep post-fit inference tied to a reconstructed, validated target.
   Optimizer moments, mini-batch gradient variance and parameter covariance
   are not interchangeable. Qualify population output and uncertainty
   separately for any new approximate model.

**Gate:** sampling/gradient oracles where unbiasedness is claimed,
full-target comparisons to deterministic fits, seed/resume tests and repeated
training studies recording failures, bias, variance, time and memory.
Unsupported dependence/sampler combinations fail explicitly.

### Phase 17 — graph accelerators

**Status:** planned, not implemented. **Prerequisite:** Phase 14 and the
specific activation/penalty/training phases being accelerated.
**Deliverable:** optional graph-likelihood acceleration with the same
qualified statistical target as CPU execution.

1. Profile before offloading. Select graph forward/derivative contractions,
   batched evaluations and eligible solves by measured total cost, including
   host/device transfers and mode/outer iterations. Existing Henderson CUDA
   kernels do not establish graph-family backend support.
2. Define operation-local device ownership, explicit device selection,
   precision and CPU/backend eligibility. Preserve signed/block family
   curvature, relationship precision order and proper latent RR semantics.
   Reuse factor/solve interfaces instead of requiring a device dense inverse.
3. Admit device memory before transfer/allocation, accounting for graph,
   derivative, solver and temporary buffers. Support bounded batching only
   where it preserves the statistical objective. Errors and requested-backend
   unavailability must be visible; any automatic CPU fallback needs an
   explicit policy and recorded reason, not silent success.
4. Establish FP64 parity first. Gate mixed precision through score residuals,
   conditional curvature, log-determinant error and refinement bounds.
   Faster arithmetic with materially different likelihood or convergence
   is not a qualified implementation.
5. Specify deterministic/reproducibility guarantees for reductions and RNG.
   Retain portable serializable fit objects; reconstruct device state for
   post-fit operations rather than saving device pointers/handles.
6. Test supported device/toolchain combinations, absence of hardware,
   insufficient memory, precision failures, cancellation and repeated
   operation cleanup. Report synchronized end-to-end times, peak device/host
   memory and crossover sizes; do not infer speedup from asynchronous launch time.

**Gate:** CPU/device objective, derivative, mode, accepted-update and output
parity for all advertised combinations, plus real hardware evidence.
Unsupported devices/platforms remain unqualified; CPU defaults are unchanged.

### Phase 18 — extended-model qualification and publication

**Status:** planned, not implemented. **Prerequisite:** all preceding
extensions included in the proposed release envelope.
**Deliverable:** an evidence-backed cross-feature support contract.

1. Cross the new axes with family/link, fixed/additive/random coordinates,
   joint/relationship/RR maps, inference target and route/backend. Test
   representative interactions, not only one successful example per phase.
2. Separate numerical derivative error, nonsmooth-region error, Laplace
   approximation, stochastic optimization error and accelerator precision.
   Use dense/root/quadrature/penalized references appropriate to each claim.
3. Run standalone recovery, coverage where justified, held-out-group
   prediction and scaling studies. Record seeds, selection rules, every
   failed fit and Monte Carlo uncertainty. Do not generalize smooth
   unpenalized coverage to selected/penalized/nonsmooth network weights.
4. Run full regressions, per-case runtime gates, installed-source/package/
   manual/vignette checks, fresh-process/checkpoint reconstruction and
   supported-platform/device CI. No mandatory `lme4`/`nlme` or accelerator
   runtime is introduced into ordinary CPU tests.
5. Publish current versus planned capability tables, model/optimizer/
   approximation labels, migration examples, executable vignettes,
   numerical options and benchmark evidence together. A dependency-complete
   feature is not automatically release-qualified.

**Dependency order:** Phase 12 establishes the baseline. Phase 13 and smooth
Phase 14 can proceed independently after it. Phase 15 uses the admitted
Phase 14 envelope; Phase 16 depends on Phase 14 and the objective/activation
features it uses. Phase 17 can follow deterministic Phase 14 without waiting
for stochastic training, but must qualify any Phase 15/16 combinations it
later advertises. Phase 18 qualifies the selected combined envelope.

**Exit gate:** each advertised combination meets its statistical, numerical,
resource and platform criteria, with explicit rejection of other requests.
No phase bypasses classical-model stability or silently changes ordinary
Gaussian/PQL defaults.

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

Begin Phase 12 by closing the outstanding platform/manual and classical
stability gates against the Phase 11 baseline. Then design and implement
Phase 13's activation contract and Phase 14's larger smooth-graph envelope,
followed by the explicitly dependent penalty, stochastic and accelerator
work. The Phase 13--18 entries are plans, not currently supported features.
Phase 7 remains skipped and unqualified until its Phase 12 closure evidence
is obtained. Preserve the public formula/method/start/control separation and
ordinary Gaussian/PQL behavior throughout.
