# Matern spatial covariance structure

Matern spatial covariance structure. The returned dimensionless
covariance shape is intended for use inside `vsm`, which supplies the
single overall variance scale.

## Usage

``` r
maternm(x, range = NULL, nu = 0.5, fixed = c(FALSE, FALSE),
        distance = NULL)
```

## Arguments

- x:

  Numeric coordinate vector or matrix/data frame whose rows give spatial
  coordinates.

- range:

  Positive starting range parameter. If NULL, the median positive
  pairwise distance is used.

- nu:

  Positive starting Matern smoothness parameter.

- fixed:

  Logical of length one or two; if length one it is recycled to range
  and nu.

- distance:

  Optional finite symmetric q by q distance matrix with zero diagonal,
  aligned to the unique spatial locations.

## Details

For distance \\d\\, the correlation is
\$\$K(d)=\frac{2^{1-\nu}}{\Gamma(\nu)}z^{\nu}K\_{\nu}(z),\qquad
z=\frac{\sqrt{2\nu}d}{\phi},\$\$ where \\\phi\\ is `range`. Both range
and smoothness are positive and are optimized on log scales. The
implementation uses the generic R evaluator and central factor-level
numerical derivatives. Repeated coordinate rows are mapped to the same
covariance level.

## Value

A list containing the incidence/design matrix in `Z` and a compiled
CovarianceFactor v2 descriptor in `covFactor`, for use inside `vsm`.

## See also

[`sar`](https://covaruber.github.io/sommer/reference/sar.md),
[`car`](https://covaruber.github.io/sommer/reference/car.md),
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md).

## Examples

``` r
if (FALSE) { # \dontrun{
vsm(maternm(cbind(xcoord,ycoord)), ism(genotype))
vsm(maternm(time), ism(id))
} # }
```
