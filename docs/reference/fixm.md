# fixed indication matrix

`fixm` creates a square matrix with 3's in the diagnals and
off-diagonals to quickly specify a fixed constraint in the Gtc argument
of the [`vsm`](https://covaruber.github.io/sommer/reference/vsm.md)
function.

## Usage

``` r
fixm(x, reps=NULL)
```

## Arguments

- x:

  integer specifying the number of traits to be fitted for a given
  random effect.

- reps:

  integer specifying the number of times the matrix should be repeated
  in a list format to provide easily the constraints in complex models
  that use the ds(), us() or cs() structures.

## Value

- \$res:

  a matrix or a list of matrices with the constraints to be provided in
  the Gtc argument of the
  [`vsm`](https://covaruber.github.io/sommer/reference/vsm.md) function.

## References

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative
traits using the R package sommer. PLoS ONE 11(6):
doi:10.1371/journal.pone.0156744

## Author

Giovanny Covarrubias-Pazaran

## Examples

``` r
fixm(4)
#>      [,1] [,2] [,3] [,4]
#> [1,]    3    3    3    3
#> [2,]    3    3    3    3
#> [3,]    3    3    3    3
#> [4,]    3    3    3    3
fixm(4,2)
#> [[1]]
#>      [,1] [,2] [,3] [,4]
#> [1,]    3    3    3    3
#> [2,]    3    3    3    3
#> [3,]    3    3    3    3
#> [4,]    3    3    3    3
#> 
#> [[2]]
#>      [,1] [,2] [,3] [,4]
#> [1,]    3    3    3    3
#> [2,]    3    3    3    3
#> [3,]    3    3    3    3
#> [4,]    3    3    3    3
#> 
```
