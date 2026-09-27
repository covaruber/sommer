# data frame to matrix

This function takes a matrix that is in data frame format and transforms
it into a matrix. Other packages that allows you to obtain an additive
relationship matrix from a pedigree is the \`pedigreemm\` package.

## Usage

``` r
dfToMatrix(x, row="Row",column="Column",
             value="Ainverse", returnInverse=FALSE, 
             bend=1e-6)
```

## Arguments

- x:

  ginv element, output from the Ainverse function.

- row:

  name of the column in x that indicates the row in the original
  relationship matrix.

- column:

  name of the column in x that indicates the column in the original
  relationship matrix.

- value:

  name of the column in x that indicates the value for a given row and
  column in the original relationship matrix.

- returnInverse:

  a TRUE/FALSE value indicating if the inverse of the x matrix should be
  computed once the data frame x is converted into a matrix.

- bend:

  a numeric value to add to the diagonal matrix in case matrix is
  singular for inversion.

## Value

- K:

  pedigree transformed in a relationship matrix.

- Kinv:

  inverse of the pedigree transformed in a relationship matrix.

## References

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative
traits using the R package sommer. PLoS ONE 11(6):
doi:10.1371/journal.pone.0156744

## Author

Giovanny Covarrubias-Pazaran

## Examples

``` r
library(Matrix)
m <- matrix(1:9,3,3)
m <- tcrossprod(m)

mdf <- as.data.frame(as.table(m))
mdf
#>   Var1 Var2 Freq
#> 1    A    A   66
#> 2    B    A   78
#> 3    C    A   90
#> 4    A    B   78
#> 5    B    B   93
#> 6    C    B  108
#> 7    A    C   90
#> 8    B    C  108
#> 9    C    C  126

dfToMatrix(mdf, row = "Var1", column = "Var2", 
            value = "Freq",returnInverse=FALSE )
#> $K
#> 3 x 3 sparse Matrix of class "dgCMatrix"
#>                
#> [1,] 66  78  90
#> [2,] 78  93 108
#> [3,] 90 108 126
#> 
#> $Kinv
#> NULL
#> 

```
