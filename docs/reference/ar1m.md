# First-order autoregressive covariance structure

First-order autoregressive covariance structure. The returned
dimensionless covariance shape is intended for use inside `vsm`, which
supplies the single overall variance scale.

## Usage

``` r
ar1m(x, rho = 0.30, fixed = FALSE,
     variance = c("homogeneous", "heterogeneous"), values = NULL)
```

## Arguments

- x:

  Ordered factor or design defining q ordered covariance levels.

- rho:

  Starting AR(1) correlation in (-1,1).

- fixed:

  For homogeneous variance, a logical indicating whether rho is fixed.
  For heterogeneous variance, a logical vector of length q: rho followed
  by q-1 variance-ratio parameters.

- variance:

  Whether the AR correlation has homogeneous or heterogeneous marginal
  variances.

- values:

  Optional q positive starting variances when
  `variance = "heterogeneous"`.

## Details

The homogeneous AR(1) shape is \$\$K\_{ij}=\rho^{\|i-j\|}.\$\$
Heterogeneous mode uses \$\$K=D R(\rho)D,\$\$ with the first variance
fixed as the relative-scale reference. The working correlation
coordinate is \\\eta=\operatorname{atanh}(\rho)\\, ensuring
\\\|\rho\|\<1\\.

## Value

A list containing the incidence/design matrix in `Z` and a compiled
CovarianceFactor v2 descriptor in `covFactor`, for use inside `vsm`.

## See also

[`ar2m`](https://covaruber.github.io/sommer/reference/ar2m.md),
[`ar3m`](https://covaruber.github.io/sommer/reference/ar3m.md),
[`toeplitzm`](https://covaruber.github.io/sommer/reference/toeplitzm.md),
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md).

## Examples

``` r
if (FALSE) { # \dontrun{
vsm(ar1m(time), ism(id))
} # }
```
