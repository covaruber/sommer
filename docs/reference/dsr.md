# Diagonal covariance structure for `mmer` and `vsr`

Creates a diagonal variance-component structure for use by the `mmer`
covariance-model interface, typically inside
[`vsr`](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md).
Each column or level represented by `x` receives its own diagonal
variance component, while all off-diagonal covariance components are
constrained to zero.

## Usage

``` r
dsr(x)
```

## Arguments

- x:

  A factor, character vector, numeric vector, or design/incidence matrix
  defining the covariance dimension. Factor and character inputs are
  expanded to incidence columns. A matrix is used directly.

## Details

`dsr()` is part of the covariance-structure interface used by
`mmer`/`vsr`; it is distinct from
[`dsm`](https://covaruber.github.io/sommer/reference/dsm.md), which
belongs to the newer `mmes`/`vsm` CovarianceFactor interface.

The function returns a design matrix \\Z\\ and a covariance-parameter
constraint matrix `thetaC`. If \\q\\ columns or levels are represented
by `x`, then \$\$\mathrm{thetaC}=I_q.\$\$ Consequently, the associated
covariance matrix has the form
\$\$\Sigma=\mathrm{diag}(\sigma_1^2,\ldots,\sigma_q^2),\$\$ so that the
\\q\\ variances are estimated separately and all covariances are fixed
to zero.

For a factor or character vector, missing observations are retained
through the construction of the incidence matrix. If only one
non-missing level is present, a one-column incidence matrix is
constructed explicitly.

For a numeric non-factor input, `x` is treated as a single design column
rather than expanded into factor levels.

## Value

A list with components:

- `Z`: the design/incidence matrix associated with the covariance
  dimension.

- `thetaC`: a diagonal matrix defining the diagonal covariance-component
  pattern used by `vsr`.

## See also

[`vsr`](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md),
[`usr`](https://covaruber.github.io/sommer/reference/usr.md),
[`csr`](https://covaruber.github.io/sommer/reference/csr.md),
[`atr`](https://covaruber.github.io/sommer/reference/atr.md),
[`dsm`](https://covaruber.github.io/sommer/reference/dsm.md),
[`mmer`](https://covaruber.github.io/sommer/reference/mmer.md)

## Examples

``` r
# Diagonal covariance structure across levels
D <- dsr(factor(c("A","B","A","C")))
D$Z
#>   A B C
#> 1 1 0 0
#> 2 0 1 0
#> 3 1 0 0
#> 4 0 0 1
#> attr(,"assign")
#> [1] 1 1 1
#> attr(,"contrasts")
#> attr(,"contrasts")$dummy
#> [1] "contr.treatment"
#> 
D$thetaC
#>   A B C
#> A 1 0 0
#> B 0 1 0
#> C 0 0 1
```
