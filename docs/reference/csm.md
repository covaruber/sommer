# Compound-symmetry covariance structure

Compound-symmetry covariance structure. The returned dimensionless
covariance shape is intended for use inside `vsm`, which supplies the
single overall variance scale.

## Usage

``` r
csm(x, rho = 0.10, fixed = FALSE,
    variance = c("homogeneous", "heterogeneous"), values = NULL)
```

## Arguments

- x:

  Variable defining q covariance levels.

- rho:

  Starting common off-diagonal correlation.

- fixed:

  For homogeneous variance, a logical indicating whether rho is fixed.
  For heterogeneous variance, a logical vector of length q: rho followed
  by q-1 variance-ratio parameters.

- variance:

  Whether the correlation has homogeneous or heterogeneous marginal
  variances.

- values:

  Optional q positive starting variances when
  `variance = "heterogeneous"`.

## Details

With `variance = "homogeneous"`, `csm` uses unit diagonal and a common
off-diagonal correlation, \$\$K\_{ii}=1, \quad K\_{ij}=\rho\\ (i\ne
j).\$\$ With `variance = "heterogeneous"`, it uses \$\$K=D C(\rho)D,\$\$
where the first variance is the reference and the remaining q-1
variances are estimated as positive ratios. Positive definiteness
requires \\-1/(q-1)\<\rho\<1\\. The product-level variance remains owned
by `vsm`.

## Value

A list containing the incidence/design matrix in `Z` and a compiled
CovarianceFactor v2 descriptor in `covFactor`, for use inside `vsm`.

## See also

[`corgm`](https://covaruber.github.io/sommer/reference/corgm.md),
[`dsm`](https://covaruber.github.io/sommer/reference/dsm.md),
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md).

## Examples

``` r
if (FALSE) { # \dontrun{
vsm(csm(environment), ism(genotype))
} # }
```
