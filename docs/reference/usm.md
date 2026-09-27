# Unstructured positive-definite covariance structure

Unstructured positive-definite covariance structure. The returned
dimensionless covariance shape is intended for use inside `vsm`, which
supplies the single overall variance scale.

## Usage

``` r
usm(x, theta = NULL, fixed = NULL)
```

## Arguments

- x:

  Variable or design defining q covariance levels.

- theta:

  Optional positive-definite q by q starting covariance matrix.

- fixed:

  Logical vector of length q(q+1)/2 - 1 controlling the estimable
  normalized-Cholesky parameters.

## Details

The covariance shape is parameterized as \$\$K=LL^{\mathsf T},\$\$ with
lower-triangular \\L\\ and \\L\_{11}=1\\. Diagonal elements after the
first are positive and represented on log scales; lower off-diagonal
elements are unrestricted. This guarantees positive definiteness without
repeatedly repairing an unconstrained covariance matrix.

## Value

A list containing the incidence/design matrix in `Z` and a compiled
CovarianceFactor v2 descriptor in `covFactor`, for use inside `vsm`.

## See also

[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md),
[`mmes`](https://covaruber.github.io/sommer/reference/mmes.md).

## Examples

``` r
if (FALSE) { # \dontrun{
vsm(usm(environment), ism(genotype))
} # }
```
