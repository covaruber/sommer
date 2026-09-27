# Native-scale covariance parameters and standard errors

Returns descriptor-defined native covariance parameters with
delta-method standard errors computed from the fitted
covariance-parameter uncertainty matrix.

## Usage

``` r
covparams_mmes_se(object, term = NULL, rel_step = 1e-6)
```

## Arguments

- object:

  A fitted object of class `"mmes"`.

- term:

  Optional covariance-term name or numeric index. By default, all random
  and residual covariance terms are returned.

- rel_step:

  Positive relative step used to numerically differentiate each
  CovarianceFactor native reporting callback.

## Details

The fitted `theta_se` matrix describes uncertainty in the reported
`covPar` coordinates. For each covariance structure, `covparams_mmes_se`
numerically evaluates the Jacobian of its descriptor's native reporting
callback and applies the multivariate delta method, \$\$V\_{native} = J
V\_{covPar} J^{\mathsf T}.\$\$ Consequently, coupled transformations use
the complete covariance information. For example, the standard error of
a [`dsm()`](https://covaruber.github.io/sommer/reference/dsm.md) level
variance accounts for uncertainty and covariance in both the overall
variance and its variance ratio. Likewise,
[`usm()`](https://covaruber.github.io/sommer/reference/usm.md) standard
errors account for all relevant normalized Cholesky parameters.

The transformation is available for every built-in CovarianceFactor
family supported by
[`covparams_mmes()`](https://covaruber.github.io/sommer/reference/covparams_mmes.md);
user-defined
[`ownm()`](https://covaruber.github.io/sommer/reference/ownm.md)
structures use their descriptor-provided `native_report` callback.

Parameters marked fixed in the covariance descriptor contribute no
independent sampling uncertainty: their rows and columns in `theta_se`
are excluded before transformation. A native quantity involving both
fixed and free parameters can still have a nonzero standard error
through its free parameters. Both fixed and estimated native quantities
are returned.

The same table is stored automatically in `object$covParNativeSE` by
[`mmes()`](https://covaruber.github.io/sommer/reference/mmes.md).

## Value

A data frame containing `term`, `factor`, `parameter`, `estimate`,
`StdError`, and `Zratio`. A zero standard error is reported for a fully
fixed native quantity, with `Zratio=NA`.

## See also

[`covparams_mmes`](https://covaruber.github.io/sommer/reference/covparams_mmes.md),
[`mmes`](https://covaruber.github.io/sommer/reference/mmes.md),
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- mmes(Yield ~ Env,
            random = ~ vsm(dsm(Env), ism(Name)),
            rcov = ~ units, data = DT_example)
covparams_mmes_se(fit)
fit$covParNativeSE
} # }
```
