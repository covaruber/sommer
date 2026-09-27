# Third-order autoregressive covariance structure

Third-order autoregressive covariance structure. The returned
dimensionless covariance shape is intended for use inside `vsm`, which
supplies the single overall variance scale.

## Usage

``` r
ar3m(x, pacf = c(0.20, 0.10, 0.05), fixed = NULL,
     variance = c("homogeneous", "heterogeneous"), values = NULL)
```

## Arguments

- x:

  Ordered factor defining covariance levels.

- pacf:

  Length-three vector of starting partial autocorrelations, each
  strictly between -1 and 1.

- fixed:

  For homogeneous variance, a logical vector of length three. For
  heterogeneous variance, a logical vector of length q+2: three PACFs
  followed by q-1 variance-ratio parameters.

- variance:

  Whether the AR correlation has homogeneous or heterogeneous marginal
  variances.

- values:

  Optional q positive starting variances when
  `variance = "heterogeneous"`.

## Details

A stationary AR(3) covariance is parameterized through three partial
autocorrelations. The PACF parameterization guarantees stationarity
while allowing unconstrained working coordinates. Heterogeneous mode
uses \$\$K=D R D\$\$ and estimates q-1 positive variance ratios in
addition to the PACFs.

## Value

A list containing the incidence/design matrix in `Z` and a compiled
CovarianceFactor v2 descriptor in `covFactor`, for use inside `vsm`.

## See also

[`ar1m`](https://covaruber.github.io/sommer/reference/ar1m.md),
[`ar2m`](https://covaruber.github.io/sommer/reference/ar2m.md),
[`toeplitzm`](https://covaruber.github.io/sommer/reference/toeplitzm.md),
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md).

## Examples

``` r
if (FALSE) { # \dontrun{
vsm(ar3m(time), ism(id))
} # }
```
