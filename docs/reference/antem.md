# Antedependence covariance structure

Antedependence covariance structure. The returned dimensionless
covariance shape is intended for use inside `vsm`, which supplies the
single overall variance scale.

## Usage

``` r
antem(x, order = 1L, beta = NULL, innovations = NULL, fixed = NULL)
```

## Arguments

- x:

  Ordered variable defining q covariance levels.

- order:

  Antedependence order, 1 \<= order \< q.

- beta:

  Optional starting regression coefficients for allowed subdiagonals.

- innovations:

  Optional q positive innovation variances.

- fixed:

  Logical vector controlling all beta coefficients and q-1
  innovation-variance ratios.

## Details

The modified-Cholesky representation is \$\$Ty=e,\qquad
\mathrm{Cov}(e)=D,\qquad K=T^{-1}DT^{-\mathsf T}.\$\$ The matrix \\T\\
is unit lower triangular, with nonzero regression coefficients only
within the requested number of preceding levels. The first innovation
variance is the scale reference; the remaining innovation variances are
positive ratios.

## Value

A list containing the incidence/design matrix in `Z` and a compiled
CovarianceFactor v2 descriptor in `covFactor`, for use inside `vsm`.

## See also

[`ar1m`](https://covaruber.github.io/sommer/reference/ar1m.md),
[`toeplitzm`](https://covaruber.github.io/sommer/reference/toeplitzm.md),
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md).

## Examples

``` r
if (FALSE) { # \dontrun{
vsm(antem(time, order=2), ism(id))
} # }
```
