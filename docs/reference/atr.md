# Selected-level diagonal covariance structure for `mmer` and `vsr`

Creates a diagonal covariance-component selector for specified columns
or levels of a covariance dimension. It is intended for the `mmer`
covariance-model interface, typically inside
[`vsr`](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md).

## Usage

``` r
atr(x, levs)
```

## Arguments

- x:

  A factor, character vector, numeric vector, or design/incidence matrix
  defining the covariance dimension. Factor and character inputs are
  expanded to incidence columns. A matrix is used directly.

- levs:

  Character vector identifying the columns or levels whose diagonal
  entries are selected. If omitted, all available column names are
  selected.

## Details

`atr()` constructs a design matrix \\Z\\ together with a diagonal
`thetaC` selector. Let \\q\\ be the number of columns of \\Z\\, and
define \$\$ a_j = \begin{cases} 1, & \text{if level }j\text{ is included
in \code{levs}},\\ 0, & \text{otherwise}. \end{cases} \$\$ The returned
constraint matrix is
\$\$\mathrm{thetaC}=\mathrm{diag}(a_1,\ldots,a_q).\$\$

Thus `atr()` is useful when a covariance component is to be associated
only with a specified subset of levels or design columns. Levels not
selected in `levs` receive zero on the corresponding diagonal of
`thetaC`.

This function belongs to the `mmer`/`vsr` covariance-structure
interface. The related
[`atm`](https://covaruber.github.io/sommer/reference/atm.md) function
belongs to the newer `mmes`/`vsm` CovarianceFactor interface and uses a
different parameterization.

For a matrix input, its column names define the available levels. For
factor or character input, the incidence-matrix column names define the
levels. If `levs` is missing, all available columns are selected.

## Value

A list with components:

- `Z`: the design/incidence representation of `x`.

- `thetaC`: a diagonal zero/one matrix selecting the levels specified by
  `levs`.

## See also

[`vsr`](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md),
[`dsr`](https://covaruber.github.io/sommer/reference/dsr.md),
[`usr`](https://covaruber.github.io/sommer/reference/usr.md),
[`csr`](https://covaruber.github.io/sommer/reference/csr.md),
[`atm`](https://covaruber.github.io/sommer/reference/atm.md),
[`mmer`](https://covaruber.github.io/sommer/reference/mmer.md)

## Examples

``` r
# Select levels B and C
A <- atr(factor(c("A","B","A","C")), levs = c("B","C"))
A$Z
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
A$thetaC
#>   A B C
#> A 0 0 0
#> B 0 1 0
#> C 0 0 1
```
