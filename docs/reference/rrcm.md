# Reduced-rank covariance structure

Reduced-rank covariance structure. The returned dimensionless covariance
shape is intended for use inside `vsm`, which supplies the single
overall variance scale.

## Usage

``` r
rrcm(x, k = 1L, loadings = NULL, fixed = NULL)
```

## Arguments

- x:

  Variable defining q covariance levels.

- k:

  Reduced rank, satisfying 1 \<= k \< q.

- loadings:

  Optional finite q by k starting loading matrix.

- fixed:

  Logical vector controlling the free loading parameters.

## Details

The reduced-rank structure uses \$\$M=\Lambda\Lambda^{\mathsf
T}+I_q,\qquad K=M/M\_{11}.\$\$ The rank-\\k\\ term captures the dominant
covariance pattern and the identity term provides a common isotropic
remainder, keeping the covariance strictly positive definite for
precision-based Henderson calculations. The leading \\k\times k\\
loading block is lower triangular for rotational identification, and its
diagonal loadings are positive. `rrcm` is compiled through the generic
CovarianceFactor callback interface rather than requiring model-specific
solver code.

## Value

A list containing the incidence/design matrix in `Z` and a compiled
CovarianceFactor v2 descriptor in `covFactor`, for use inside `vsm`.

## See also

[`fam`](https://covaruber.github.io/sommer/reference/fam.md),
[`ownm`](https://covaruber.github.io/sommer/reference/ownm.md),
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md).

## Examples

``` r
if (FALSE) { # \dontrun{
vsm(rrcm(environment, k=2), ism(genotype))
} # }
```
