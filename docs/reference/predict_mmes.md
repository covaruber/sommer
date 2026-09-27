# Predict form of a LMM fitted with mmes

`predict` method for class `"mmes"`.

## Usage

``` r
# S3 method for class 'mmes'
predict(object, Dtable=NULL, D, ...)
```

## Arguments

- object:

  a mixed model of class `"mmes"`

- Dtable:

  a table specifying the terms to be included or averaged.

  An "include" term means that the model matrices for that fixed or
  random effect is filled with 1's for the positions where column names
  and row names match.

  An "include and average" term means that the model matrices for that
  fixed or random effect is filled with 1/#1's in that row.

  An "average" term alone means that all rows for such fixed or random
  effect will be filled with 1/#levels in the effect.

  If a term is not considered "include" or "average" is then totally
  ignored in the BLUP and SE calculation.

  The default rule to invoke when the user doesn't provide the Dtable is
  to include and average all terms that match the argument D.

- D:

  a character string specifying the variable used to extract levels for
  the rows of the D matrix and its construction. Alternatively, the D
  matrix (of class dgCMatrix) specifying the matrix to be used for the
  predictions directly.

- ...:

  Further arguments to be passed.

## Details

This function allows to produce predictions specifying those variables
that define the margins of the hypertable to be predicted (argument D).
Predictions are obtained for each combination of values of the specified
variables that is present in the data set used to fit the model. See
vignettes for more details.

Standard errors are exact for any `computeCi` setting (0, 1, or 2) and
do not require `object$Ci`: `Var(D %*% bu) = D C.inv() D.t()` is instead
obtained by solving `C x = d` for each row `d` of `D` against the stored
coefficient matrix `object$C`, which
[`mmes()`](https://covaruber.github.io/sommer/reference/mmes.md) always
returns regardless of `computeCi`. This is cheaper than forming the
complete inverse (`computeCi=2`) and, unlike the Takahashi
selected-inverse subset (`computeCi=1`), is exact for arbitrary linear
combinations, not just diagonal PEVs.

For predicted values the pertinent design matrices X and Z together with
BLUEs (b) and BLUPs (u) are multiplied and added together.

predicted.value equal Xb + Zu.1 + ... + Zu.n

For computing standard errors for predictions the parts of the
coefficient matrix:

C11 equal (X.t() V.inv() X).inv()

C12 equal 0 - \[(X.t() V.inv() X).inv() X.t() V.inv() G Z\]

C22 equal PEV equal G - \[Z.t() G\[V.inv() - (V.inv() X X.t() V.inv() X
V.inv() X)\]G Z.t()\]

In practive C equals ( W.t() V.inv() W ).inv()

when both fixed and random effects are present in the inclusion set. If
only fixed and random effects are included, only the respective terms
from the SE for fixed or random effects are calculated.

## Value

- pvals:

  the table of predictions according to the specified arguments.

- vcov:

  the variance covariance for the predictions.

- D:

  the model matrix for predictions as defined in Welham et al.(2004).

- Dtable:

  the table specifying the terms to include and terms to be averaged.

## References

Welham, S., Cullis, B., Gogel, B., Gilmour, A., and Thompson, R. (2004).
Prediction in linear mixed models. Australian and New Zealand Journal of
Statistics, 46, 325 - 347.

## Author

Giovanny Covarrubias-Pazaran

## See also

[`predict`](https://rdrr.io/r/stats/predict.html),
[`mmes`](https://covaruber.github.io/sommer/reference/mmes.md)

## Examples

``` r
data(DT_yatesoats, package="enhancer")
DT <- DT_yatesoats
m3 <- mmes(fixed=Y ~ V + N + V:N ,
           random = ~ B + B:MP,
           rcov=~units,
           data = DT)
#> Solver selected: ldlt
#> OpenMP available: up to 8 threads.
#> OpenMP active: parallel SLQ log-determinant probes (8 probes).
#>   Working-coordinate trust scaling applied: alpha=0.830121
#> iteration    LogLik     wall    cpu(sec)   restrained   EM weight      pivot
#>     1      -13.8039   21:18:6      0           0      1      4.81656
#>     2      -11.5674   21:18:6      0           0      0.813615      6.80191
#>     3      -11.4975   21:18:6      0           0      0.661969      6.76374
#>     4      -11.4965   21:18:6      0           0      0.538588      6.80222
#>     5      -11.4965   21:18:6      0           0      0.438203      6.807

m3 <- postPEV(m3, mode = 2)
 
#############################
## predict means for nitrogen
#############################
Dt <- m3$Dtable; Dt
#>     type           term include average
#> 1  fixed              1   FALSE   FALSE
#> 2  fixed              V   FALSE   FALSE
#> 3  fixed              N   FALSE   FALSE
#> 4  fixed            V:N   FALSE   FALSE
#> 5 random    vsm(ism(B))   FALSE   FALSE
#> 6 random vsm(ism(B:MP))   FALSE   FALSE
# first fixed effect just average
Dt[1,"average"] = TRUE
# second fixed effect include
Dt[2,"include"] = TRUE
# third fixed effect include and average
Dt[3,"include"] = TRUE
Dt[3,"average"] = TRUE
Dt
#>     type           term include average
#> 1  fixed              1   FALSE    TRUE
#> 2  fixed              V    TRUE   FALSE
#> 3  fixed              N    TRUE    TRUE
#> 4  fixed            V:N   FALSE   FALSE
#> 5 random    vsm(ism(B))   FALSE   FALSE
#> 6 random vsm(ism(B:MP))   FALSE   FALSE

pp=predict(object=m3, Dtable=Dt, D="N")
pp$pvals
#>       N predicted.value std.error
#> 0     0        78.16667  13.31593
#> 0.2 0.2        82.79167  13.99152
#> 0.4 0.4        86.83333  13.99152
#> 0.6 0.6        89.37500  13.99152

#############################
## predict means for variety
#############################

Dt <- m3$Dtable; Dt
#>     type           term include average
#> 1  fixed              1   FALSE   FALSE
#> 2  fixed              V   FALSE   FALSE
#> 3  fixed              N   FALSE   FALSE
#> 4  fixed            V:N   FALSE   FALSE
#> 5 random    vsm(ism(B))   FALSE   FALSE
#> 6 random vsm(ism(B:MP))   FALSE   FALSE
# first fixed effect include
Dt[1,"include"] = TRUE
# second fixed effect just average
Dt[2,"average"] = TRUE
# third fixed effect include and average
Dt[3,"include"] = TRUE
Dt[3,"average"] = TRUE
Dt
#>     type           term include average
#> 1  fixed              1    TRUE   FALSE
#> 2  fixed              V   FALSE    TRUE
#> 3  fixed              N    TRUE    TRUE
#> 4  fixed            V:N   FALSE   FALSE
#> 5 random    vsm(ism(B))   FALSE   FALSE
#> 6 random vsm(ism(B:MP))   FALSE   FALSE

pp=predict(object=m3, Dtable=Dt, D="V")
pp$pvals
#>                     V predicted.value std.error
#> GoldenRain GoldenRain        103.8889  7.671521
#> Marvellous Marvellous        103.8889  7.671521
#> ictory         ictory        103.8889  7.671521

#############################
## predict means for nitrogen:variety
#############################
# prediction matrix D based on (equivalent to classify in asreml)
Dt <- m3$Dtable; Dt
#>     type           term include average
#> 1  fixed              1   FALSE   FALSE
#> 2  fixed              V   FALSE   FALSE
#> 3  fixed              N   FALSE   FALSE
#> 4  fixed            V:N   FALSE   FALSE
#> 5 random    vsm(ism(B))   FALSE   FALSE
#> 6 random vsm(ism(B:MP))   FALSE   FALSE
# first fixed effect include and average
Dt[1,"include"] = TRUE
Dt[1,"average"] = TRUE
# second fixed effect include and average
Dt[2,"include"] = TRUE
Dt[2,"average"] = TRUE
# third fixed effect include and average
Dt[3,"include"] = TRUE
Dt[3,"average"] = TRUE
Dt
#>     type           term include average
#> 1  fixed              1    TRUE    TRUE
#> 2  fixed              V    TRUE    TRUE
#> 3  fixed              N    TRUE    TRUE
#> 4  fixed            V:N   FALSE   FALSE
#> 5 random    vsm(ism(B))   FALSE   FALSE
#> 6 random vsm(ism(B:MP))   FALSE   FALSE

pp=predict(object=m3, Dtable=Dt, D="N:V")
pp$pvals
#>                               N:V predicted.value std.error
#> N0:VGoldenRain     N0:VGoldenRain        80.00000  9.106759
#> N0.2:VGoldenRain N0.2:VGoldenRain        84.62500  8.477256
#> N0.4:VGoldenRain N0.4:VGoldenRain        88.66667  8.477256
#> N0.6:VGoldenRain N0.6:VGoldenRain        91.20833  8.477256
#> N0:VMarvellous     N0:VMarvellous        82.22222  7.871438
#> N0.2:VMarvellous N0.2:VMarvellous        86.84722  7.470608
#> N0.4:VMarvellous N0.4:VMarvellous        90.88889  7.470608
#> N0.6:VMarvellous N0.6:VMarvellous        93.43056  7.470608
#> N0:VVictory           N0:VVictory        77.16667  7.871438
#> N0.2:VVictory       N0.2:VVictory        81.79167  7.470608
#> N0.4:VVictory       N0.4:VVictory        85.83333  7.470608
#> N0.6:VVictory       N0.6:VVictory        88.37500  7.470608
```
