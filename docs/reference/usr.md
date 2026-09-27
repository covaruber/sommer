# Unstructured covariance structure for `mmer` and `vsr`

Creates an unstructured variance-covariance component pattern for use by
the `mmer` covariance-model interface, typically inside
[`vsr`](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md).

## Usage

``` r
usr(x)
```

## Arguments

- x:

  A factor, character vector, numeric vector, or design/incidence matrix
  defining the covariance dimension. Factor and character inputs are
  expanded to incidence columns. A matrix is used directly.

## Details

`usr()` is part of the covariance-structure interface used by
`mmer`/`vsr`; it is distinct from
[`usm`](https://covaruber.github.io/sommer/reference/usm.md), which
belongs to the newer `mmes`/`vsm` CovarianceFactor interface.

After constructing the design matrix \\Z\\, the function calls
`unsm(q)`, where \\q\\ is the number of columns of \\Z\\, to create the
`thetaC` pattern for an unstructured covariance matrix. Thus the
associated covariance model permits a separate variance for every level
and a separate covariance for every pair of levels: \$\$ \Sigma =
\begin{bmatrix} \sigma_1^2 & \sigma\_{12} & \cdots & \sigma\_{1q}\\
\sigma\_{12} & \sigma_2^2 & \cdots & \sigma\_{2q}\\ \vdots & \vdots &
\ddots & \vdots\\ \sigma\_{1q} & \sigma\_{2q} & \cdots & \sigma_q^2
\end{bmatrix}. \$\$ The exact component coding is stored in the matrix
returned by `unsm(q)` and passed to the `vsr` machinery through
`thetaC`.

For a factor or character vector, missing observations are retained
through the construction of the incidence matrix. For a numeric
non-factor input, `x` is treated as a single design column.

## Value

A list with components:

- `Z`: the design/incidence matrix associated with the covariance
  dimension.

- `thetaC`: the unstructured covariance-component pattern produced by
  [`unsm()`](https://covaruber.github.io/sommer/reference/unsm.md).

## See also

[`vsr`](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md),
[`dsr`](https://covaruber.github.io/sommer/reference/dsr.md),
[`csr`](https://covaruber.github.io/sommer/reference/csr.md),
[`atr`](https://covaruber.github.io/sommer/reference/atr.md),
[`usm`](https://covaruber.github.io/sommer/reference/usm.md),
[`mmer`](https://covaruber.github.io/sommer/reference/mmer.md)

## Examples

``` r
# Unstructured covariance pattern across three levels
U <- usr(factor(c("A","B","A","C")))
U$Z
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
U$thetaC
#>   A B C
#> A 1 2 2
#> B 2 1 2
#> C 2 2 1
```
