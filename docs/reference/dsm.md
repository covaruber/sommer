# Diagonal heterogeneous covariance structure

Diagonal heterogeneous covariance structure. The returned dimensionless
covariance shape is intended for use inside `vsm`, which supplies the
single overall variance scale.

## Usage

``` r
dsm(x, values = NULL, fixed = NULL, theta = NULL)
```

## Arguments

- x:

  Variable or design defining the covariance levels.

- values:

  Optional positive starting variances, one per level.

- fixed:

  Logical vector of length q-1 indicating which variance-ratio
  parameters are fixed.

- theta:

  Optional q by q diagonal starting covariance matrix; if supplied, its
  diagonal replaces values.

## Details

For \\q\\ levels, `dsm` uses the first variance as the internal scale
reference and estimates \\q-1\\ positive variance ratios,
\$\$K=\mathrm{diag}(1,r_2,\ldots,r_q).\$\$ The ratios are optimized on
log scales and reported on their natural positive scale.

## Value

A list containing the incidence/design matrix in `Z` and a compiled
CovarianceFactor v2 descriptor in `covFactor`, for use inside `vsm`.

## See also

[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md),
[`mmes`](https://covaruber.github.io/sommer/reference/mmes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
vsm(dsm(environment), ism(genotype))
} # }
```
