# Nonlinear mixed-model effects: phased implementation plan

## Goal

Add a native, likelihood-based nonlinear predictor for mixed models. The
predictor must support learned nonlinear fixed effects and, in a later phase,
nonlinear random effects. A one-hidden-layer network such as

```
eta_i = relu(x_i W + c) v + d + z_i u
```

is an initial target, not a special-case implementation. Other supported
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

## Source review: current fitting architecture

The existing code has useful seams, but none currently represents a nonlinear
predictor in estimated parameters.

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

## Proposed model and API contract

Use one declarative predictor graph for fixed and random nonlinear components.
Do not settle public syntax until a small API prototype proves that formulas,
parameter names, dimensions, prediction, and serialization work together.

The API should express, at minimum:

- Input columns or existing prepared terms.
- An ordered list of layers/functions, each with a declared output size and
  activation or mean function.
- Which parameter blocks are fixed and which are group-specific random
  parameters.
- Parameter names, initialization, and supported constraints.
- Optional conventional sommer random effects added to the nonlinear
  predictor.

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

**Deliverable:** a private, small prototype of a declarative one-hidden-layer
predictor and an agreed public model equation/API before touching Gaussian
solvers.

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
2. Add a native predictor evaluator for a bounded initial function set, such
   as identity, ReLU, tanh, and logistic. Implement the initial network
   architecture from composition of these primitives, not hard-coded
   one-off equations.
3. Optimize nonlinear fixed parameters jointly with covariance and residual
   parameters against the Gaussian ML objective. Keep the existing random
   predictor linear in this phase, so its established conditional mode and
   Hessian remain valid.
4. Use a numerical outer gradient initially only if needed as a reference or
   fallback. It must be explicit, bounded by validated parameter/workspace
   limits, and reported. Do not claim it scales to large networks.
5. Store the graph, parameter layout, training preprocessing, and fitted
   parameters in the result. Add new-data prediction by evaluating the same
   graph, not by manufacturing `model.matrix()` columns for learned layers.
6. Route the request through an explicit model-path planner. Unsupported
   controls must error before fitting; do not silently select a linear fit.

**Initial restrictions:** Gaussian ML, independent diagonal residuals,
   Henderson/native Laplace fitting, and conventional additive Gaussian
   random terms. Reject REML, PQL, direct fitting, structured residuals,
   unsupported family mixtures, and unsupported control combinations until
   each has a proved implementation.

**Tests and gate:**

- Identity/linear graph matches existing Gaussian ML estimates, likelihood,
  fitted values, and covariance estimates within numerical tolerance.
- Fixed ReLU graph matches a direct least-squares or independent optimizer
  reference on a small full-rank dataset.
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
4. Keep categorical block likelihoods, zero-inflated models, family mixtures,
   and Gaussian structured residual blocks out until their cross-derivative
   and indexing contracts are tested.

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

### Phase 5 — scalability and release qualification

**Deliverable:** a qualified feature envelope, not a blanket claim of general
neural-network support.

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

**Stop condition:** stop expanding when the documented supported envelope
passes numerical, regression, prediction, memory, and package qualification
gates. Add new activations or model structures only with their own derivative
and recovery tests.

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

## Recommended first implementation slice

Implement Phase 0, then Phase 1 for one-hidden-layer Gaussian ML with a fixed
ReLU predictor and conventional linear sommer random effects. Keep this slice
Henderson-only, use explicit ML, preserve ordinary fixed/random formulas for
their current terms, and add an opt-in nonlinear specification. This proves
the parameter graph, optimizer, prediction, and diagnostics without modifying
the established Gaussian AI engines or prematurely changing random-mode
curvature.
