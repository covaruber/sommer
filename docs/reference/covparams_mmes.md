# Extract native-scale covariance parameters from an mmes fit

Returns covariance parameters in the native, model-facing scale defined
by each CovarianceFactor descriptor, rather than the unconstrained
working coordinates or normalized variance ratios used internally by the
optimizer.

## Usage

``` r
covparams_mmes(object, term = NULL)
```

## Arguments

- object:

  A fitted object of class `"mmes"`.

- term:

  Optional covariance-term name or numeric index. By default, all random
  and residual covariance terms are returned.

## Details

Every CovarianceFactor descriptor carries a native reporting callback.
Built-in structures use it to report meaningful quantities such as
level-specific variances for
[`dsm()`](https://covaruber.github.io/sommer/reference/dsm.md),
variances and covariances for
[`usm()`](https://covaruber.github.io/sommer/reference/usm.md), and
variance-correlation parameters for correlation structures. Structural
zeros and internal coordinates such as log standard deviations,
normalized Cholesky coefficients, and variance ratios are not returned
where the covariance model provides a more direct native representation.

All built-in CovarianceFactor families provide native reporters,
including identity, diagonal, unstructured, AR(1–3), compound-symmetry,
moving-average, general-correlation, factor-analytic, antedependence,
reduced-rank, Matern, Toeplitz, SAR, and CAR structures. Heterogeneous
AR and compound-symmetry models report their dependence parameter
together with level-specific variances rather than an overall variance
plus variance ratios.

Both fixed and estimated covariance parameters are included. For an
arbitrary Kronecker product, the product-level variance scale is
absorbed by the first covariance-shaping factor. Subsequent factors are
necessarily reported on their normalized relative scale because the
product has only one identifiable overall scale.

User-defined
[`ownm()`](https://covaruber.github.io/sommer/reference/ownm.md)
structures can supply a `native_report` callback. If none is supplied,
their overall variance and naturally transformed parameter values are
returned.

Use
[`covparams_mmes_se()`](https://covaruber.github.io/sommer/reference/covparams_mmes_se.md)
when native-scale delta-method standard errors are also required.

## Value

A data frame with columns `term`, `factor`, `parameter`, and `estimate`.

## See also

[`covparams_mmes_se`](https://covaruber.github.io/sommer/reference/covparams_mmes_se.md),
[`mmes`](https://covaruber.github.io/sommer/reference/mmes.md),
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md),
[`ownm`](https://covaruber.github.io/sommer/reference/ownm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- mmes(Yield ~ Env,
            random = ~ vsm(dsm(Env), ism(Name)),
            rcov = ~ units, data = DT_example)
covparams_mmes(fit)
fit$covParNative
} # }
```
