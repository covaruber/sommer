# Classical nonlinear mixed models

The nonlinear model interface supports parametric covariate effects, expressed
as bounded formula expressions with declared `parameter()` models or as native
exponential, logistic, and saturation curves. Curves are declared with
`curves=list(name=nonlinear_curve(...))` and referenced by `curve_term("name")`.

The existing Gaussian linear, linear mixed, generalized linear, generalized
linear mixed, and nonlinear mixed model engines remain independent of this
interface. Nonlinear mixed models support additive Gaussian random effects and
the admitted random curve-parameter covariance structures. Generalized
nonlinear models use the supported family likelihoods and integrations.

Nonlinear fits use maximum likelihood and explicit numerical controls.
Prediction may target conditional, zero-deviation, or population means where
the selected model supports them. Uncertainty, conditional covariance, and
prediction intervals are explicit post-fit operations.

## Release qualification

Before release, run the retained classical nonlinear and ordinary model
regressions, clean native builds, package examples, vignettes, and package
checks. Compare retained likelihoods and predictions with independent
references. Cross-platform and PDF toolchain qualification remain separate
release gates and are not implied by local test results.
