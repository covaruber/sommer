# Factor-analytic covariance structure

Factor-analytic covariance structure. The returned dimensionless
covariance shape is intended for use inside `vsm`, which supplies the
single overall variance scale.

## Usage

``` r
fam(x, k = 1L, loadings = NULL, specific = NULL, fixed = NULL)
```

## Arguments

- x:

  Variable defining q covariance levels.

- k:

  Factor-analytic rank, satisfying 1 \<= k \< q.

- loadings:

  Optional finite q by k starting loading matrix.

- specific:

  Optional q positive starting specific variances.

- fixed:

  Logical vector controlling all loading and specific-variance-ratio
  parameters.

## Details

The factor-analytic shape is based on \$\$M=\Lambda\Lambda^{\mathsf
T}+\Psi,\$\$ where \\\Psi\\ is diagonal positive, followed by
\\K=M/M\_{11}\\. For rotational identification, the leading \\k\times
k\\ loading block is lower triangular and its diagonal loadings are
positive. The first specific variance is an internal reference;
remaining specific variances are positive ratios. This removes the
otherwise redundant internal scale because `vsm` supplies the overall
\\\sigma^2\\.

## Value

A list containing the incidence/design matrix in `Z` and a compiled
CovarianceFactor v2 descriptor in `covFactor`, for use inside `vsm`.

## See also

[`rrcm`](https://covaruber.github.io/sommer/reference/rrcm.md),
[`usm`](https://covaruber.github.io/sommer/reference/usm.md),
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md).

## Examples

``` r
if (FALSE) { # \dontrun{
vsm(fam(environment, k=2), ism(genotype))
} # }
```
