# Identity covariance structure

Identity covariance structure. The returned dimensionless covariance
shape is intended for use inside `vsm`, which supplies the single
overall variance scale.

## Usage

``` r
ism(x)
```

## Arguments

- x:

  Factor, character vector, numeric design, or matrix defining
  covariance levels.

## Details

`ism` returns \\K=I\\. It has no covariance-shape parameters; the only
unknown covariance scale is the \\\sigma^2\\ supplied by `vsm`. It is
commonly used as the final term of `vsm` to indicate an independent main
effect.

## Value

A list containing the incidence/design matrix in `Z` and a compiled
CovarianceFactor v2 descriptor in `covFactor`, for use inside `vsm`.

## See also

[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md),
[`mmes`](https://covaruber.github.io/sommer/reference/mmes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
vsm(ism(id))
vsm(dsm(environment), ism(genotype))
} # }
```
