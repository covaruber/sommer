# User-defined covariance structure

User-defined covariance structure. The returned dimensionless covariance
shape is intended for use inside `vsm`, which supplies the single
overall variance scale.

## Usage

``` r
ownm(x, K = NULL, fun = NULL, par = numeric(), fixed = NULL,
  dfun = NULL, par_names = NULL, native_report = NULL)
```

## Arguments

- x:

  Variable or design defining q covariance levels.

- K:

  Optional known positive-definite q by q covariance matrix. If
  supplied, the structure is fixed.

- fun:

  Optional function fun(par) returning a finite q by q covariance
  matrix.

- par:

  Starting parameter vector for fun.

- fixed:

  Logical vector of length par indicating fixed parameters.

- dfun:

  Optional derivative function dfun(par, k) returning the derivative of
  the raw covariance matrix with respect to parameter k.

- par_names:

  Optional names for the parameters.

- native_report:

  Optional function with arguments `scale`, `par`, `factor`, and
  `absorb_scale`. It must return a named numeric vector giving the
  fitted parameters in the model's preferred native scale. This callback
  is used by
  [`covparams_mmes()`](https://covaruber.github.io/sommer/reference/covparams_mmes.md).

## Details

There are two modes. With `K`, `ownm` creates a fixed covariance factor.
With `fun`, the matrix is evaluated at every optimizer point. The
returned matrix is symmetrized and normalized by its first diagonal
element so that `vsm` retains the unique overall variance scale. If
`dfun` is supplied, its derivative is normalized analytically using the
quotient rule; otherwise a central numerical derivative of the
normalized factor is used. This is the general extension mechanism for
covariance structures that do not require a native C++ evaluator. When
`native_report` is omitted, the overall variance and naturally
transformed `par` values are reported.

## Value

A list containing the incidence/design matrix in `Z` and a compiled
CovarianceFactor v2 descriptor in `covFactor`, for use inside `vsm`.

## See also

[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md),
[`covparams_mmes`](https://covaruber.github.io/sommer/reference/covparams_mmes.md),
[`rrcm`](https://covaruber.github.io/sommer/reference/rrcm.md),
[`maternm`](https://covaruber.github.io/sommer/reference/maternm.md),
[`sar`](https://covaruber.github.io/sommer/reference/sar.md),
[`car`](https://covaruber.github.io/sommer/reference/car.md).

## Examples

``` r
if (FALSE) { # \dontrun{
K <- matrix(c(1,.3,.3,1),2,2)
vsm(ownm(group, K=K), ism(id))

cfun <- function(p) {
  r <- tanh(p[1]); matrix(c(1,r,r,1),2,2)
}
native <- function(scale, par, factor, absorb_scale=TRUE) {
  c(variance=scale, correlation=tanh(par[1]))
}
vsm(ownm(group, fun=cfun, par=0, native_report=native), ism(id))
} # }
```
