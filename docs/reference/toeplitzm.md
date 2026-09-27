# General positive-definite Toeplitz correlation structure

General positive-definite Toeplitz correlation structure. The returned
dimensionless covariance shape is intended for use inside `vsm`, which
supplies the single overall variance scale.

## Usage

``` r
toeplitzm(x, pacf = NULL, fixed = NULL)
```

## Arguments

- x:

  Ordered factor or design defining q covariance levels.

- pacf:

  Optional q-1 starting reflection coefficients / partial
  autocorrelations. Defaults to 0.10.

- fixed:

  Logical vector of length q-1.

## Details

A full q by q Toeplitz correlation matrix is parameterized by q-1
reflection coefficients. Each PACF is mapped from an unconstrained
working coordinate with `tanh`; the reflection coefficients are
converted to a stable autoregressive representation and then to the
implied Toeplitz correlation sequence. This parameterization preserves
positive definiteness while allowing every lag correlation to vary
indirectly.

## Value

A list containing the incidence/design matrix in `Z` and a compiled
CovarianceFactor v2 descriptor in `covFactor`, for use inside `vsm`.

## See also

[`ar1m`](https://covaruber.github.io/sommer/reference/ar1m.md),
[`ar2m`](https://covaruber.github.io/sommer/reference/ar2m.md),
[`ar3m`](https://covaruber.github.io/sommer/reference/ar3m.md),
[`mam`](https://covaruber.github.io/sommer/reference/mam.md),
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md).

## Examples

``` r
if (FALSE) { # \dontrun{
vsm(toeplitzm(time), ism(id))
} # }
```
