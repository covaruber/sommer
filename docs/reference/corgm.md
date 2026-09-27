# General positive-definite correlation structure

General positive-definite correlation structure. The returned
dimensionless covariance shape is intended for use inside `vsm`, which
supplies the single overall variance scale.

## Usage

``` r
corgm(x, theta = NULL, fixed = NULL)
```

## Arguments

- x:

  Variable defining q covariance levels.

- theta:

  Optional positive-definite q by q starting covariance or correlation
  matrix.

- fixed:

  Logical vector of length q(q-1)/2.

## Details

An unrestricted SPD correlation matrix is represented using a
unit-diagonal lower factor \\A\\. The intermediate matrix
\\S=AA^{\mathsf T}\\ is standardized to correlation scale,
\$\$K\_{ij}=S\_{ij}/\sqrt{S\_{ii}S\_{jj}}.\$\$ The unrestricted working
parameters are the strict lower-triangular entries of \\A\\.

## Value

A list containing the incidence/design matrix in `Z` and a compiled
CovarianceFactor v2 descriptor in `covFactor`, for use inside `vsm`.

## See also

[`csm`](https://covaruber.github.io/sommer/reference/csm.md),
[`usm`](https://covaruber.github.io/sommer/reference/usm.md),
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md).

## Examples

``` r
if (FALSE) { # \dontrun{
vsm(corgm(environment), ism(genotype))
} # }
```
