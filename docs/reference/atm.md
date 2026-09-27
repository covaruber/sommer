# Selected-level diagonal covariance structure

Selected-level diagonal covariance structure. The returned dimensionless
covariance shape is intended for use inside `vsm`, which supplies the
single overall variance scale.

## Usage

``` r
atm(x, levs, values = NULL, fixed = NULL)
```

## Arguments

- x:

  Variable defining all available levels.

- levs:

  Character names or numeric positions of the levels retained in this
  covariance factor.

- values:

  Optional positive starting variances for the selected levels.

- fixed:

  Logical vector of length length(levs)-1 controlling the selected
  variance ratios.

## Details

`atm` is a selected-level version of `dsm`. The design is restricted to
`levs`; observations outside those levels receive zeros for this factor.
The first selected level is the scale reference and the remaining
selected levels are positive variance ratios. This representation is
most naturally used for random effects; residual factors must still
assign exactly one covariance-product coordinate to every observation.

## Value

A list containing the incidence/design matrix in `Z` and a compiled
CovarianceFactor v2 descriptor in `covFactor`, for use inside `vsm`.

## See also

[`dsm`](https://covaruber.github.io/sommer/reference/dsm.md),
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md),
[`mmes`](https://covaruber.github.io/sommer/reference/mmes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
vsm(atm(environment, c("E1","E3")), ism(genotype))
} # }
```
