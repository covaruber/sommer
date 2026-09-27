# User-specified covariance-component structure for `mmer` and `vsr`

Combines a design/incidence representation with a user-supplied
covariance-component constraint matrix for use by the `mmer`
covariance-model interface, typically inside
[`vsr`](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md).

## Usage

``` r
csr(x, mm)
```

## Arguments

- x:

  A factor, character vector, numeric vector, or design/incidence matrix
  defining the covariance dimension. Factor and character inputs are
  expanded to incidence columns. A matrix is used directly.

- mm:

  A square covariance-component constraint matrix. Its dimensions must
  correspond to the columns of the design matrix generated from `x`. The
  matrix is returned as `thetaC`, with its row and column names set to
  the column names of the resulting design matrix.

## Details

`csr()` provides the general constraint-matrix interface underlying
custom covariance patterns for `mmer`/`vsr`. Unlike
[`dsr`](https://covaruber.github.io/sommer/reference/dsr.md) and
[`usr`](https://covaruber.github.io/sommer/reference/usr.md), it does
not construct a predefined diagonal or unstructured pattern. Instead,
the supplied matrix `mm` determines the covariance-component structure
through `thetaC`.

Conceptually, `thetaC` identifies which elements of the covariance
matrix are associated with estimated variance-covariance components and
which constraints are shared among elements. The interpretation of the
entries therefore follows the `thetaC` convention used by
[`vsr`](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md).

When `x` is a matrix it is used directly as \\Z\\. Otherwise, factors
and character vectors are converted to incidence matrices, while numeric
non-factor inputs are treated as a single design column.

The function assigns the column names of \\Z\\ to both dimensions of
`mm`; therefore `mm` should be conformable with the number of columns of
the resulting design matrix.

## Value

A list with components:

- `Z`: the design/incidence matrix associated with the covariance
  dimension.

- `thetaC`: the user-supplied covariance-component constraint matrix,
  with dimensions named according to `Z`.

## See also

[`vsr`](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md),
[`dsr`](https://covaruber.github.io/sommer/reference/dsr.md),
[`usr`](https://covaruber.github.io/sommer/reference/usr.md),
[`atr`](https://covaruber.github.io/sommer/reference/atr.md),
[`mmer`](https://covaruber.github.io/sommer/reference/mmer.md)

## Examples

``` r
# A custom two-level covariance-component pattern
x <- factor(c("A","B","A","B"))
mm <- matrix(c(1, 2,
               2, 1), 2, 2, byrow = TRUE)
C <- csr(x, mm)
C$Z
#>   A B
#> 1 1 0
#> 2 0 1
#> 3 1 0
#> 4 0 1
#> attr(,"assign")
#> [1] 1 1
#> attr(,"contrasts")
#> attr(,"contrasts")$dummy
#> [1] "contr.treatment"
#> 
C$thetaC
#>   A B
#> A 1 2
#> B 2 1
```
