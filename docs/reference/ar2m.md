# Second-order autoregressive covariance structure

Second-order autoregressive covariance structure. The returned
dimensionless covariance shape is intended for use inside `vsm`, which
supplies the single overall variance scale.

## Usage

``` r
ar2m(x, pacf = c(0.20, 0.10), fixed = NULL,
     variance = c("homogeneous", "heterogeneous"), values = NULL)
```

## Arguments

- x:

  Ordered factor defining covariance levels.

- pacf:

  Length-two vector of starting partial autocorrelations, each strictly
  between -1 and 1.

- fixed:

  For homogeneous variance, a logical vector of length two. For
  heterogeneous variance, a logical vector of length q+1: two PACFs
  followed by q-1 variance-ratio parameters.

- variance:

  Whether the AR correlation has homogeneous or heterogeneous marginal
  variances.

- values:

  Optional q positive starting variances when
  `variance = "heterogeneous"`.

## Details

A stationary AR(2) covariance is parameterized through reflection
coefficients / partial autocorrelations. Each PACF is mapped from an
unconstrained working coordinate by `tanh`, then converted to stable AR
coefficients by the Levinson/Durbin recursion. Heterogeneous mode uses
\$\$K=D R D\$\$ and estimates q-1 positive variance ratios in addition
to the PACFs.

## Value

A list containing the incidence/design matrix in `Z` and a compiled
CovarianceFactor v2 descriptor in `covFactor`, for use inside `vsm`.

## See also

[`ar1m`](https://covaruber.github.io/sommer/reference/ar1m.md),
[`ar3m`](https://covaruber.github.io/sommer/reference/ar3m.md),
[`toeplitzm`](https://covaruber.github.io/sommer/reference/toeplitzm.md),
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md).

## Examples

``` r
if (FALSE) { # \dontrun{
vsm(ar2m(time), ism(id))
} # }
```
