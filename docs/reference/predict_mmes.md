# Predict form of a LMM fitted with mmes

`predict` method for class `"mmes"`.

## Usage

``` r
# S3 method for class 'mmes'
predict(object, Dtable=NULL, D, levels=NULL, sed=FALSE,
  pairwise=FALSE, adjust="none", df=Inf, PEV=TRUE, VarU=TRUE, ...)
```

## Arguments

- object:

  a mixed model of class `"mmes"`

- Dtable:

  a table specifying the terms to be included or averaged.

  An "include" term uses cells matching each prediction row. Categorical
  random effects are mapped within each fitted term, including fitted
  cells without records. Numeric fixed terms retain their evaluated values.

  An "include and average" term averages the included cells. For fixed
  categorical terms the count includes reference cells absorbed in the intercept.

  An "average" term alone averages over selected factor levels. Numeric
  fixed terms are evaluated at their covariate settings before averaging.

  If a term is not considered "include" or "average" is then totally
  ignored in the BLUP and SE calculation.

  The default rule to invoke when the user doesn't provide the Dtable is
  to include and average all terms that match the argument D.

  The `levels` column is a list-column. `NULL` means all factor levels and
  default numeric settings. The returned table records resolved numeric
  settings and can be reused.

- D:

  a character string specifying the variable used to extract levels for
  the rows of the D matrix and its construction. Alternatively, the D
  matrix (of class dgCMatrix) specifying the matrix to be used for the
  predictions directly.

- levels:

  An optional named list whose names match `Dtable$term`. Fixed-factor
  interactions use data cell labels, such as `"Victory:0.2"`. Structured
  random terms use fitted cell labels, such as `"CA.2011:genotype1"`.
  For an included categorical classify term, levels select prediction rows
  in the supplied order and can request fitted levels present only in `Gu`.
  For other terms, they restrict inclusion or averaging.

  A term with one numeric variable accepts a scalar, for example
  `levels=list("Env:cov"=5)`. Numeric or mixed terms also accept a named
  per-variable list, for example
  `levels=list("Env:cov"=list(Env=c("CA.2011", "CA.2012"), cov=5))`.
  Numeric settings must be finite scalar values. Use the exact term names
  in the fitted Dtable. Restrictions are term-specific: restricting `V`
  does not automatically restrict `V:N`.

- sed:

  Return standard errors of differences and their summary when `TRUE`.

- pairwise:

  `TRUE` for all pairwise differences, or a single reference prediction level.

- adjust:

  Multiplicity adjustment passed to `stats::p.adjust`.

- df:

  Degrees of freedom: `Inf`, a positive number, `"residual"`,
  `"satterthwaite"`, or `"kr"`. Small-sample methods apply only to
  fixed-effect contrasts; contrasts involving random effects remain normal-based.

- PEV:

  Logical, default `TRUE`. Return prediction-error covariance for the
  requested linear combinations as `PEV` and its diagonal as `pvals$pev`.
  When `FALSE`, omit these additional outputs; `vcov` and `std.error`
  remain available for compatibility.

- VarU:

  Logical, default `TRUE`. Return sampling covariance of the random
  contribution to each requested linear combination as `VarU`, its diagonal
  as `pvals$var.u`, and whole-predictor sampling covariance as `sampling.vcov`.
  Set `FALSE` to skip the extra relationship solves. Requires a Gaussian
  Henderson fit; use `VarU=FALSE` for PQL working-model predictions.

- ...:

  Further arguments to be passed.

## Details

This function allows to produce predictions specifying those variables
that define the margins of the hypertable to be predicted (argument D).
Predictions are obtained for each combination of values of the specified
variables that is present in the data set used to fit the model. See
vignettes for more details.

Categorical structured random terms use their fitted covariance-coordinate
and main-effect levels to construct genetic weights, even for cells without
observations. Unsupported random designs that cannot be mapped for all rows
produce a warning; use an explicit `D` matrix for these designs.

Numeric fixed terms, including factor-by-covariate interactions, use the
fitted formula and contrasts. Numeric variables default to their means in
the retained fitting data, unless overridden. Factor levels are equally
weighted when averaged: an averaged `Env:cov` term with one slope per
environment uses `mean(cov)/nlevels(Env)` for each slope.

Write $D=[L_b,L_u]$. The target is $t=L_b\beta+L_u u$ and its predictor is
$\hat t=L_b\hat\beta+L_u\hat u$. All prediction covariance outputs refer to
the final returned `D`, after `Dtable` and `levels` have been applied. They
have one row and column per requested prediction, not per model effect.

Let $C_0$ be the original-scale Henderson coefficient matrix,
$Q=(C_0^{-1})_{uu}$, and $B=(C_0^{-1})_{bb}$. Then

$$
\mathrm{PEV}=\operatorname{Var}(t-\hat t)=D C_0^{-1}D',
$$
$$
\mathrm{VarU}=L_u(G-Q)L_u',
$$
$$
\mathrm{sampling.vcov}=L_b B L_b'+L_u(G-Q)L_u'.
$$

The existing `vcov` equals PEV and `std.error` is its diagonal square root.
The fixed BLUE and random BLUP have zero sampling cross-covariance, but
their estimation/prediction errors generally do not. Subtracting full-target
PEV from $L_uGL_u'$ is therefore incorrect when fixed effects are included.
Fixed-only predictions have zero random `VarU`, but generally nonzero
sampling covariance and PEV. For random-only predictions,
$L_uGL_u'=\mathrm{PEV}+\mathrm{VarU}$.

Computations batch full, random-only and fixed-only contrasts through one
Henderson factorization, solving against their transposes and applying
`Cscale`. They do not use `Ci` or construct/invert observation covariance
$V$. Prior covariance products use $G_k=\Sigma_k\otimes A_k$ and sparse
relationship-precision solves. New fits retain precisions in original
level order; older ordinary fits attempt input reconstruction. Rotated
and factor-score effects are mapped to public coordinates. Direct and
matrix-free `solveOnly` fits lack the required stored coefficient system.

These are exact known-parameter covariances and plug-in approximations at
estimated variance components, omitting parameter-estimation uncertainty.
Additional matrices require quadratic storage in the prediction count.
SEDs and pairwise tests retain their PEV-based meaning.

## Value

- pvals:

  the table of predictions according to the specified arguments.

- vcov:

  Full target prediction-error covariance, always returned.

- PEV:

  Identical to `vcov` when requested; its diagonal is `pvals$pev`.

- VarU:

  Sampling covariance of the requested random contribution;
  its diagonal is `pvals$var.u`.

- sampling.vcov:

  Sampling covariance of the whole requested predictor, when `VarU=TRUE`.

- D:

  the model matrix for predictions as defined in Welham et al.(2004).

- Dtable:

  the table specifying included and averaged terms and resolved levels/settings.

- sed, avsed, pairwise:

  Requested standard errors of differences, their summary, and pairwise comparisons.

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

The `levels` argument and the `Dtable$levels` list-column are equivalent.
This example selects nitrogen levels in the requested order:

```r
Dt <- m3$Dtable
Dt$average[Dt$term %in% c("1", "V", "V:N")] <- TRUE
Dt$include[Dt$term %in% c("N", "V:N")] <- TRUE

p <- predict(m3, D="N", Dtable=Dt,
             levels=list(N=c("0.6", "0.2")))
p$pvals
p$Dtable$levels[[match("N", p$Dtable$term)]]

Dt$levels[[match("N", Dt$term)]] <- c("0.6", "0.2")
p2 <- predict(m3, D="N", Dtable=Dt)
all.equal(p$D, p2$D)
#> [1] TRUE
```

The covariance matrices are projected through the requested linear combination:

```r
p <- predict(m3, D="N", Dtable=Dt, PEV=TRUE, VarU=TRUE)
p$PEV
p$VarU
p$sampling.vcov
randomColumns <- seq.int(nrow(m3$b)+1L, nrow(m3$bu))
randomD <- p$D[,randomColumns,drop=FALSE]
full <- postVarU(m3, mode=2)
all.equal(unname(p$VarU),
          unname(as.matrix(randomD %*% full$VarU %*% t(randomD))))
```
