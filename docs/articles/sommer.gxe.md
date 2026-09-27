# Fitting genotype by environment models in sommer

The sommer package was developed to provide R users with a flexible
univariate and multivariate linear mixed-model solver. The current
package provides two complementary REML formulations. The
[`mmes()`](https://covaruber.github.io/sommer/reference/mmes.md)
function uses Henderson’s mixed model equations and a sparse Average
Information REML algorithm: the mixed-model coefficient matrix is
factorized with sparse LDLT methods, selected inverse elements required
by REML calculations can be obtained with Takahashi sparse-inverse
recursions, and the factorization is reused for mixed-model solutions
and variance-parameter derivatives. The
[`mmer()`](https://covaruber.github.io/sommer/reference/mmer.md)
function uses the marginal-covariance MNR formulation, in which REML
calculations are organized around the observation covariance matrix
$`V`$ and the REML projection matrix $`P`$. The relative efficiency of
the two formulations depends on model dimensions, sparsity, and
covariance structure rather than only on whether $`p>n`$ or $`n>p`$. The
core numerical algorithms are coded in C++ using Armadillo and Eigen.
The package allows users to specify flexible variance-covariance
structures for random and residual effects and to obtain quantities such
as REML variance-covariance estimates, BLUPs, BLUEs, residuals, fitted
values, and prediction error variance information.

The purpose of this vignette is to show how to fit different genotype by
environment (GxE) models using the sommer package:

1.  Single environment model
2.  Multienvironment model: Main effect model
3.  Multienvironment model: Diagonal model (DG)
4.  Multienvironment model: Compund symmetry model (CS)
5.  Multienvironment model: Unstructured model (US)
6.  Multienvironment model: Random regression model (RR)
7.  Multienvironment model: Other covariance structures for GxE
8.  Multienvironment model: Finlay-Wilkinson regression
9.  Multienvironment model: Factor analytic (reduced rank) model (FA)
10. Two stage analysis

When the breeder decides to run a trial and apply selection in a single
environment (whether because the amount of seed is a limitation or
there’s no availability for a location) the breeder takes the risk of
selecting material for a target population of environments (TPEs) using
an environment that is not representative of the larger TPE. Therefore,
many breeding programs try to base their selection decision on
multi-environment trial (MET) data. Models could be adjusted by adding
additional information like spatial information, experimental design
information, etc. In this tutorial we will focus mainly on the
covariance structures for GxE and the incorporation of relationship
matrices for the genotype effect.

## 1) Single environment model

A single-environment model is the one that is fitted when the breeding
program can only afford one location, leaving out the possible
information available from other environments. This will be used to
further expand to GxE models.

``` r
library(sommer)
```

    ## Loading required package: Matrix

    ## Loading required package: MASS

    ## Loading required package: crayon

    ## Loading required package: enhancer

``` r
data(DT_example, package="enhancer")
DT <- DT_example
A <- A_example

Ai <- solve(A)
Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Ai, "inverse")=TRUE

ansSingle <- mmes(Yield~1,
              random= ~ vsm(ism(Name), Gu=Ai),
              rcov= ~ units,
              data=DT, verbose = FALSE)
```

    ## Solver selected: ldlt

``` r
summary(ansSingle)
```

    ## ============================================================
    ##          Multivariate Linear Mixed Model fit by  REML         
    ## **********************  sommer 4.4  ********************** 
    ## ============================================================
    ##          logLik      AIC      BIC Method Converge
    ## Value -81.41893 164.8379 168.0582     AI     TRUE
    ## ============================================================
    ## Variance-Covariance components:
    ##                      term factor parameter estimate StdError Zratio
    ## 1 vsm(ism(Name), Gu = Ai) sigma2    sigma2    6.533    1.532  4.263
    ## 2         vsm(ism(units)) sigma2    sigma2   13.865    1.143 12.133
    ## ============================================================
    ## Fixed effects:
    ##           Estimate Std.Error t.value
    ## Intercept    11.74        NA      NA
    ## ============================================================
    ## Use the '$' sign to access results and parameters

In this model, the random term is the germplasm effect (here called
`Name`). For the sake of example, a relationship structure among the
levels of `Name` is supplied through `Gu`. In the code above `Ai` is the
sparse precision matrix corresponding to the relationship matrix `A`,
and `attr(Ai, "inverse")=TRUE` tells sommer that the supplied matrix is
already in precision form. More generally, a non-diagonal relationship
structure can be used to model covariance among genotype levels.

## 2) MET: main effect model

A multi-environment model is the one that is fitted when the breeding
program can afford more than one location. The main effect model assumes
that GxE doesn’t exist and that the main genotype effect plus the fixed
effect for environment is enough to predict the genotype effect in all
locations of interest.

``` r
ansMain <- mmes(Yield~Env,
              random= ~ vsm(ism(Name), Gu=Ai),
              rcov= ~ units,
              data=DT, verbose = FALSE)
```

    ## Solver selected: ldlt

``` r
summary(ansMain)
```

    ## ============================================================
    ##          Multivariate Linear Mixed Model fit by  REML         
    ## **********************  sommer 4.4  ********************** 
    ## ============================================================
    ##          logLik      AIC      BIC Method Converge
    ## Value -38.73391 83.46782 93.12888     AI     TRUE
    ## ============================================================
    ## Variance-Covariance components:
    ##                      term factor parameter estimate StdError Zratio
    ## 1 vsm(ism(Name), Gu = Ai) sigma2    sigma2    4.857    1.054  4.609
    ## 2         vsm(ism(units)) sigma2    sigma2    8.107    0.673 12.047
    ## ============================================================
    ## Fixed effects:
    ##            Estimate Std.Error t.value
    ## Intercept    16.385        NA      NA
    ## EnvCA.2012   -5.688        NA      NA
    ## EnvCA.2013   -6.218        NA      NA
    ## ============================================================
    ## Use the '$' sign to access results and parameters

## 3) MET: diagonal model (DG)

A multi-environment model is fitted when observations are available from
more than one environment. The diagonal GxE model allows the genetic
variance to differ among environments while constraining the genetic
covariance between different environments to zero. In the current
[`vsm()`](https://covaruber.github.io/sommer/reference/vsm.md)
parameterization, `dsm(Env)` is a normalized diagonal covariance factor:
one product-level variance scale is estimated by
[`vsm()`](https://covaruber.github.io/sommer/reference/vsm.md) and the
remaining diagonal parameters are relative variance ratios. Thus the
model has an environment-specific genetic variance without estimating
cross-environment genetic covariances. The fixed environment effect and
the environment-specific genotype BLUP are then combined to predict
performance in each environment.

``` r
ansDG <- mmes(Yield~Env,
              random= ~ vsm(dsm(Env),ism(Name), Gu=Ai),
              rcov= ~ units,
              data=DT, verbose = FALSE)
```

    ## Solver selected: ldlt

``` r
summary(ansDG)
```

    ## ============================================================
    ##          Multivariate Linear Mixed Model fit by  REML         
    ## **********************  sommer 4.4  ********************** 
    ## ============================================================
    ##          logLik      AIC      BIC Method Converge
    ## Value -27.18128 60.36255 70.02362     AI     TRUE
    ## ============================================================
    ## Variance-Covariance components:
    ##                                term factor         parameter estimate StdError
    ## 1 vsm(dsm(Env), ism(Name), Gu = Ai)   diag variance[CA.2011]   17.468   3.7715
    ## 2 vsm(dsm(Env), ism(Name), Gu = Ai)   diag variance[CA.2012]    5.339   1.2664
    ## 3 vsm(dsm(Env), ism(Name), Gu = Ai)   diag variance[CA.2013]    7.886   1.8053
    ## 4                   vsm(ism(units)) sigma2            sigma2    4.381   0.4592
    ##   Zratio
    ## 1  4.632
    ## 2  4.216
    ## 3  4.368
    ## 4  9.540
    ## ============================================================
    ## Fixed effects:
    ##            Estimate Std.Error t.value
    ## Intercept    16.621        NA      NA
    ## EnvCA.2012   -5.958        NA      NA
    ## EnvCA.2013   -6.662        NA      NA
    ## ============================================================
    ## Use the '$' sign to access results and parameters

## 4) MET: compund symmetry model (CS)

A multi-environment model is fitted when observations are available from
more than one environment. In the code below, compound-symmetry-like GxE
behavior is represented by two independent random terms: a genotype main
effect shared across environments and an environment-by-genotype
deviation. Their covariance contribution is the sum of the two
random-effect covariance matrices. Under the identity structures used
here this implies a common covariance between observations on the same
genotype in different environments, while the GxE deviation adds
environment-specific variance. The fixed environment effect, genotype
main-effect BLUP, and genotype-by-environment BLUP together describe
performance in each environment.

``` r
E <- diag(length(unique(DT$Env)));rownames(E) <- colnames(E) <- unique(DT$Env)
Ei <- solve(E)
Ai <- solve(A)
EAi <- kronecker(Ei,Ai, make.dimnames = TRUE)
Ei <- as(as(as( Ei,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
Ai <- as(as(as( Ai,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
EAi <- as(as(as( EAi,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(Ai, "inverse")=TRUE
attr(EAi, "inverse")=TRUE
ansCS <- mmes(Yield~Env,
              random= ~ vsm(ism(Name), Gu=Ai) + vsm(ism(Env:Name), Gu=EAi),
              rcov= ~ units, 
              data=DT, verbose = FALSE)
```

    ## Solver selected: ldlt

``` r
summary(ansCS)
```

    ## ============================================================
    ##          Multivariate Linear Mixed Model fit by  REML         
    ## **********************  sommer 4.4  ********************** 
    ## ============================================================
    ##          logLik      AIC      BIC Method Converge
    ## Value -26.28507 58.57015 68.23121     AI     TRUE
    ## ============================================================
    ## Variance-Covariance components:
    ##                           term factor parameter estimate StdError Zratio
    ## 1      vsm(ism(Name), Gu = Ai) sigma2    sigma2    3.683   1.1193  3.291
    ## 2 vsm(ism(Env:Name), Gu = EAi) sigma2    sigma2    5.171   1.0038  5.151
    ## 3              vsm(ism(units)) sigma2    sigma2    4.367   0.4528  9.644
    ## ============================================================
    ## Fixed effects:
    ##            Estimate Std.Error t.value
    ## Intercept    16.496        NA      NA
    ## EnvCA.2012   -5.777        NA      NA
    ## EnvCA.2013   -6.380        NA      NA
    ## ============================================================
    ## Use the '$' sign to access results and parameters

## 5) MET: unstructured model (US)

The unstructured GxE model allows a distinct genetic variance for each
environment and a distinct genetic covariance for every pair of
environments. In the current
[`vsm()`](https://covaruber.github.io/sommer/reference/vsm.md)
architecture, `usm(Env)` is represented by a normalized Cholesky
covariance factor and
[`vsm()`](https://covaruber.github.io/sommer/reference/vsm.md) supplies
the single product-level variance scale. This parameterization
guarantees positive definiteness for admissible working parameters while
retaining the flexibility of an unstructured covariance matrix. Because
the number of covariance-shape parameters grows rapidly with the number
of environments, these models can still be statistically weakly
identified or computationally demanding. The fixed environment effect
and the correlated environment-specific genotype BLUPs are used to
predict performance in each environment.

``` r
ansUS <- mmes(Yield~Env,
              random= ~ vsm(usm(Env),ism(Name), Gu=Ai),
              rcov= ~ units,
              data=DT, verbose = FALSE)
```

    ## Solver selected: ldlt

``` r
summary(ansUS)
```

    ## ============================================================
    ##          Multivariate Linear Mixed Model fit by  REML         
    ## **********************  sommer 4.4  ********************** 
    ## ============================================================
    ##          logLik      AIC      BIC Method Converge
    ## Value -20.34919 46.69838 56.35945     AI     TRUE
    ## ============================================================
    ## Variance-Covariance components:
    ##                                term factor                   parameter estimate
    ## 1 vsm(usm(Env), ism(Name), Gu = Ai)     us           variance[CA.2011]  15.9861
    ## 2 vsm(usm(Env), ism(Name), Gu = Ai)     us covariance[CA.2012,CA.2011]   6.1695
    ## 3 vsm(usm(Env), ism(Name), Gu = Ai)     us covariance[CA.2013,CA.2011]   6.3640
    ## 4 vsm(usm(Env), ism(Name), Gu = Ai)     us           variance[CA.2012]   5.2744
    ## 5 vsm(usm(Env), ism(Name), Gu = Ai)     us covariance[CA.2013,CA.2012]   0.3745
    ## 6 vsm(usm(Env), ism(Name), Gu = Ai)     us           variance[CA.2013]   7.6896
    ## 7                   vsm(ism(units)) sigma2                      sigma2   4.3859
    ##   StdError Zratio
    ## 1    3.354 4.7669
    ## 2    1.656 3.7259
    ## 3    1.957 3.2514
    ## 4    1.266 4.1672
    ## 5    1.096 0.3417
    ## 6    1.779 4.3223
    ## 7    0.456 9.6180
    ## ============================================================
    ## Fixed effects:
    ##            Estimate Std.Error t.value
    ## Intercept    16.341        NA      NA
    ## EnvCA.2012   -5.696        NA      NA
    ## EnvCA.2013   -6.286        NA      NA
    ## ============================================================
    ## Use the '$' sign to access results and parameters

## 6) MET: random regression model

A random regression model represents environmental response with a
continuous covariate or a set of basis functions, here Legendre
polynomials of the numeric environment index. Genotype-specific
coefficients on these basis functions describe reaction norms across
environments. The covariance structure assigned to the basis
coefficients determines how the random intercept, slopes, and
higher-order coefficients vary and covary. Consequently, the number of
covariance parameters depends both on the polynomial order and on the
covariance structure placed on the basis coefficients.

``` r
library(orthopolynom)
DT$EnvN <- as.numeric(as.factor(DT$Env))

ansRR <- mmes(Yield~Env,
              random= ~ vsm(dsm(leg(EnvN,1)),ism(Name)),
              rcov= ~ units,
              data=DT, verbose = FALSE)
```

    ## Solver selected: ldlt

``` r
summary(ansRR)
```

    ## ============================================================
    ##          Multivariate Linear Mixed Model fit by  REML         
    ## **********************  sommer 4.4  ********************** 
    ## ============================================================
    ##          logLik      AIC      BIC Method Converge
    ## Value -33.84286 73.68572 83.34678     AI     TRUE
    ## ============================================================
    ## Variance-Covariance components:
    ##                                term factor      parameter estimate StdError
    ## 1 vsm(dsm(leg(EnvN, 1)), ism(Name))   diag variance[leg0]   10.392   2.0585
    ## 2 vsm(dsm(leg(EnvN, 1)), ism(Name))   diag variance[leg1]    2.080   0.6994
    ## 3                   vsm(ism(units)) sigma2         sigma2    6.296   0.5976
    ##   Zratio
    ## 1  5.048
    ## 2  2.975
    ## 3 10.535
    ## ============================================================
    ## Fixed effects:
    ##            Estimate Std.Error t.value
    ## Intercept    16.541        NA      NA
    ## EnvCA.2012   -5.832        NA      NA
    ## EnvCA.2013   -6.472        NA      NA
    ## ============================================================
    ## Use the '$' sign to access results and parameters

In addition, an unstructured, diagonal or other variance-covariance
structure can be put on top of the polynomial model:

``` r
library(orthopolynom)
DT$EnvN <- as.numeric(as.factor(DT$Env))

ansRR <- mmes(Yield~Env,
              random= ~ vsm(usm(leg(EnvN,1)),ism(Name)),
              rcov= ~ units,
              data=DT, verbose = FALSE)
```

    ## Solver selected: ldlt

``` r
summary(ansRR)
```

    ## ============================================================
    ##          Multivariate Linear Mixed Model fit by  REML         
    ## **********************  sommer 4.4  ********************** 
    ## ============================================================
    ##          logLik     AIC      BIC Method Converge
    ## Value -31.70935 69.4187 79.07977     AI     TRUE
    ## ============================================================
    ## Variance-Covariance components:
    ##                                term factor             parameter estimate
    ## 1 vsm(usm(leg(EnvN, 1)), ism(Name))     us        variance[leg0]   10.791
    ## 2 vsm(usm(leg(EnvN, 1)), ism(Name))     us covariance[leg1,leg0]   -2.428
    ## 3 vsm(usm(leg(EnvN, 1)), ism(Name))     us        variance[leg1]    2.288
    ## 4                   vsm(ism(units)) sigma2                sigma2    6.259
    ##   StdError Zratio
    ## 1   2.1617  4.992
    ## 2   0.9523 -2.550
    ## 3   0.7609  3.007
    ## 4   0.5951 10.517
    ## ============================================================
    ## Fixed effects:
    ##            Estimate Std.Error t.value
    ## Intercept    16.501        NA      NA
    ## EnvCA.2012   -5.791        NA      NA
    ## EnvCA.2013   -6.476        NA      NA
    ## ============================================================
    ## Use the '$' sign to access results and parameters

## 7) Other GxE covariance structures

Many structured covariance models can be used for GxE when their
assumptions are appropriate for the ordering or relationship among
environments. In the example below `csm(Env)` specifies a homogeneous
correlation structure for the environment covariance factor. Other
available structures, such as autoregressive models for meaningfully
ordered environments, can be substituted within the same
[`vsm()`](https://covaruber.github.io/sommer/reference/vsm.md) product
architecture.

``` r
ansAR1 <- mmes(Yield~Env,
              random= ~ vsm(csm(Env),ism(Name)),
              rcov= ~ units,
              data=DT, verbose = FALSE)
```

    ## Solver selected: ldlt

``` r
summary(ansAR1)
```

    ## ============================================================
    ##          Multivariate Linear Mixed Model fit by  REML         
    ## **********************  sommer 4.4  ********************** 
    ## ============================================================
    ##         logLik      AIC      BIC Method Converge
    ## Value -26.2851 58.57019 68.23126     AI     TRUE
    ## ============================================================
    ## Variance-Covariance components:
    ##                       term factor parameter estimate StdError Zratio
    ## 1 vsm(csm(Env), ism(Name))    csm  variance   8.8454   1.2286  7.200
    ## 2 vsm(csm(Env), ism(Name))    csm       rho   0.4156   0.1013  4.104
    ## 3          vsm(ism(units)) sigma2    sigma2   4.3697   0.4512  9.684
    ## ============================================================
    ## Fixed effects:
    ##            Estimate Std.Error t.value
    ## Intercept    16.496        NA      NA
    ## EnvCA.2012   -5.777        NA      NA
    ## EnvCA.2013   -6.381        NA      NA
    ## ============================================================
    ## Use the '$' sign to access results and parameters

## 8) Finlay-Wilkinson regression

``` r
data(DT_h2, package="enhancer")
DT <- DT_h2

## build the environmental index
ei <- aggregate(y~Env, data=DT,FUN=mean)
colnames(ei)[2] <- "envIndex"
ei$envIndex <- ei$envIndex - mean(ei$envIndex,na.rm=TRUE) # center the envIndex to have clean VCs
ei <- ei[with(ei, order(envIndex)), ]

## add the environmental index to the original dataset
DT2 <- merge(DT,ei, by="Env")

# numeric by factor variables like envIndex:Name can't be used in the random part like this
# they need to come with the vsm() structure
DT2 <- DT2[with(DT2, order(Name)), ]
mix2 <- mmes(y~ envIndex, 
             random=~ Name + vsm(ism(envIndex),ism(Name)), data=DT2,
             rcov=~vsm(dsm(Name),ism(units)),
             nIters = 50, verbose = FALSE
)
```

    ## Solver selected: ldlt

``` r
# summary(mix2)$varcomp

b=mix2$uList$`vsm(ism(envIndex), ism(Name` # adaptability (b) or genotype slopes
mu=mix2$uList$`vsm(ism(Name`# general adaptation (mu) or main effect
e=sqrt(summary(mix2)$varcomp[-c(1:2),"estimate"]) # error variance for each individual

## general adaptation (main effect) vs adaptability (response to better environments)
plot(mu[,1]~b[,1], ylab="general adaptation", xlab="adaptability")
text(y=mu[,1],x=b[,1], labels = rownames(mu), cex=0.5, pos = 1)
```

![](sommer.gxe_files/figure-html/unnamed-chunk-9-1.png)

``` r
## prediction across environments
Dt <- mix2$Dtable
Dt[1,"average"]=TRUE
Dt[2,"include"]=TRUE
Dt[3,"include"]=TRUE

mix2 <- postPEV(mix2, mode=2)
pp <- predict(mix2,Dtable = Dt, D="Name")
preds <- pp$pvals
# preds[with(preds, order(-predicted.value)), ]
## performance vs stability (deviation from regression line)
plot(preds[,2]~e, ylab="performance", xlab="stability")
text(y=preds[,2],x=e, labels = rownames(mu), cex=0.5, pos = 1)
```

![](sommer.gxe_files/figure-html/unnamed-chunk-9-2.png)

## 9) Factor analytic (reduced rank) model

When the number of environments is large, a fully unstructured genetic
covariance requires $`q(q+1)/2`$ covariance parameters for $`q`$
environments and can become weakly identified relative to the available
information. Reduced-rank and factor-analytic representations provide a
more parsimonious alternative. In the first implementation below,
[`rrm()`](https://rdrr.io/pkg/enhancer/man/rrm.html) constructs a
reduced environmental basis from `H0`;
[`usm()`](https://covaruber.github.io/sommer/reference/usm.md) estimates
the covariance among the retained factor scores, and the additional
`dsm(Env)` term models environment-specific diagonal variation. The
commented alternative uses the newer
[`rrcm()`](https://covaruber.github.io/sommer/reference/rrcm.md)
CovarianceFactor constructor directly. These representations reduce the
dimension of the cross-environment covariance model while retaining
major covariance patterns.

``` r
data(DT_h2, package="enhancer")
DT <- DT_h2
DT=DT[with(DT, order(Env)), ]
head(DT)
```

    ##          Name     Env Loc Year     Block  y
    ## 67   MSL007-B CA.2011  CA 2011 CA.2011.2  5
    ## 105  MSL007-B CA.2011  CA 2011 CA.2011.1  6
    ## 308  MSK061-4 CA.2011  CA 2011 CA.2011.2  9
    ## 393  MSK061-4 CA.2011  CA 2011 CA.2011.1 10
    ## 469 MSR169-8Y CA.2011  CA 2011 CA.2011.1 11
    ## 471     NY148 CA.2011  CA 2011 CA.2011.1 11

``` r
indNames <- na.omit(unique(DT$Name))
A <- diag(length(indNames))
rownames(A) <- colnames(A) <- indNames

# fit diagonal model first to produce H matrix
ansDG <- mmes(y~Env, henderson=TRUE,
              random=~ vsm(dsm(Env), ism(Name)),
              rcov=~units, nIters = 100,
              data=DT, verbose = FALSE)
```

    ## Solver selected: ldlt

``` r
H0 <- ansDG$uList$`vsm(dsm(Env), ism(Name))` # GxE table

# reduced rank model
ansFA <- mmes(y~Env, henderson=TRUE,
              random=~vsm( usm(rrm(Env, H = H0, nPC = 3)) , ism(Name)) + # rr
                vsm(dsm(Env), ism(Name)), # diag
              rcov=~units,
              # we recommend giving more iterations to these models
              nIters = 100, verbose = FALSE,
              # we recommend giving more EM iterations at the beggining
              data=DT)
```

    ## Solver selected: ldlt

``` r
vcFA <- ansFA$theta[[1]]
vcDG <- ansFA$theta[[2]]

loadings=with(DT, rrm(Env, nPC = 3, H = H0, returnGamma = TRUE) )$Gamma
scores <- ansFA$uList[[1]]

vcUS <- loadings %*% vcFA %*% t(loadings)
G <- vcUS + vcDG
# colfunc <- colorRampPalette(c("steelblue4","springgreen","yellow"))
# hv <- heatmap(cov2cor(G), col = colfunc(100), symm = TRUE)

uFA <- scores %*% t(loadings)
uDG <- ansFA$uList[[2]]
u <- uFA + uDG

# option 2
# ansFA2 <- mmes(y~Env, henderson=TRUE,
#               random=~vsm( rrcm(Env, 3) , ism(Name)) + # rr
#                 vsm(dsm(Env), ism(Name)), # diag
#               rcov=~units,
#               # we recommend giving more iterations to these models
#               nIters = 100, verbose = FALSE,
#               # we recommend giving more EM iterations at the beggining
#               data=DT)
# u2 <- ansFA2$uList$`vsm(rrcm(Env, 3), ism(Name` + ansFA2$uList$`vsm(dsm(Env), ism(Name`
```

For the [`rrm()`](https://rdrr.io/pkg/enhancer/man/rrm.html)
representation above, genotype BLUPs on the original environment scale
are recovered by multiplying the retained factor scores by the transpose
of the loading matrix (`Gamma`) and then adding the diagonal GxE
contribution. Likewise, the reduced-rank covariance contribution is
reconstructed as $`\Gamma\,\Sigma_f\,\Gamma^\prime`$, where $`\Sigma_f`$
is the estimated covariance among factor scores. This provides a
parsimonious approximation to a general cross-environment covariance
structure.

## 10) Two stage analysis

In two-stage analyses, a first-stage mixed model is commonly fitted to
account for experimental-design and field variation while treating the
entry effects of interest as fixed, producing adjusted entry means
(BLUEs or EMMs). Their estimation-error covariance should be carried
into the second stage rather than treating the adjusted means as
independent observations with equal precision. In the example below,
`computeCi = 2` requests the complete coefficient-matrix inverse after
fitting each first-stage model; the fixed-effect block is used to
construct the corresponding precision contribution for the second-stage
weight matrix `W`. The second-stage mixed model then analyzes the
adjusted means while retaining this first-stage precision information.

``` r
##########
## stage 1
## use mmes for dense field trials
##########
data(DT_h2, package="enhancer")
DT <- DT_h2
head(DT)
```

    ##                 Name     Env Loc Year     Block y
    ## 1            W8822-3 FL.2012  FL 2012 FL.2012.1 2
    ## 2            W8867-7 FL.2012  FL 2012 FL.2012.2 2
    ## 3           MSL007-B MO.2011  MO 2011 MO.2011.1 3
    ## 4         CO00270-7W FL.2012  FL 2012 FL.2012.2 3
    ## 5 Manistee(MSL292-A) FL.2013  FL 2013 FL.2013.2 3
    ## 6           MSM246-B FL.2012  FL 2012 FL.2012.2 3

``` r
envs <- unique(DT$Env)
BLUEL <- list()
XtXL <- list()
for(i in 1:length(envs)){
  ans1 <- mmes(y~Name-1,
               random=~Block,
               verbose=FALSE,
               computeCi = 2,
               data=droplevels(DT[which(DT$Env == envs[i]),])
  )
  ans1$Beta$Env <- envs[i]
  
  BLUEL[[i]] <- data.frame( Effect=factor(rownames(ans1$b)), 
                            Estimate=ans1$b[,1], 
                            Env=factor(envs[i]))
  # to be comparable to 1/(se^2) = 1/PEV = 1/Ci = 1/[(X'X)inv]
  XtXL[[i]] <- solve(ans1$Ci[1:nrow(ans1$b),1:nrow(ans1$b)]) 
}
```

    ## Solver selected: cholmod
    ## Solver selected: cholmod
    ## Solver selected: cholmod
    ## Solver selected: cholmod
    ## Solver selected: cholmod
    ## Solver selected: cholmod
    ## Solver selected: cholmod
    ## Solver selected: cholmod
    ## Solver selected: cholmod
    ## Solver selected: cholmod
    ## Solver selected: cholmod
    ## Solver selected: cholmod
    ## Solver selected: cholmod
    ## Solver selected: cholmod
    ## Solver selected: cholmod

``` r
DT2 <- do.call(rbind, BLUEL)
OM <- Reduce(adiag1,lapply(XtXL,as.matrix))

##########
## stage 2
## use mmes for sparse equation
##########
m <- matrix(1/var(DT2$Estimate, na.rm = TRUE))

ans2 <- mmes(Estimate~Env, henderson=TRUE,
             random=~ Effect + Env:Effect, 
             rcov = ~ vsm(
               ism(units),
               sigma2 = 1,
               fixedSigma2 = TRUE
             ),
             W=OM, 
             verbose=FALSE,
             data=DT2
)
```

    ## Solver selected: ldlt

    ## Using the weights matrix

``` r
summary(ans2)$varcomp
```

    ##                                              term factor parameter estimate
    ## 1                                vsm(ism(Effect)) sigma2    sigma2 2.076896
    ## 2                            vsm(ism(Env:Effect)) sigma2    sigma2 3.337145
    ## 3 vsm(ism(units), sigma2 = 1, fixedSigma2 = TRUE) sigma2    sigma2 1.000000
    ##    StdError    Zratio
    ## 1 0.4069953  5.102998
    ## 2 0.3167026 10.537158
    ## 3 0.0000000        NA

## Literature

Covarrubias-Pazaran G. 2016. Genome assisted prediction of quantitative
traits using the R package sommer. PLoS ONE 11(6):1-15.

Covarrubias-Pazaran G. 2018. Software update: Moving the R package
sommer to multivariate mixed models for genome-assisted prediction. doi:
<https://doi.org/10.1101/354639>

Bernardo Rex. 2010. Breeding for quantitative traits in plants. Second
edition. Stemma Press. 390 pp.

Gilmour et al. 1995. Average Information REML: An efficient algorithm
for variance parameter estimation in linear mixed models. Biometrics
51(4):1440-1450.

Henderson C.R. 1975. Best Linear Unbiased Estimation and Prediction
under a Selection Model. Biometrics vol. 31(2):423-447.

Kang et al. 2008. Efficient control of population structure in model
organism association mapping. Genetics 178:1709-1723.

Lee, D.-J., Durban, M., and Eilers, P.H.C. (2013). Efficient
two-dimensional smoothing with P-spline ANOVA mixed models and nested
bases. Computational Statistics and Data Analysis, 61, 22 - 37.

Lee et al. 2015. MTG2: An efficient algorithm for multivariate linear
mixed model analysis based on genomic information. Cold Spring Harbor.
doi: <http://dx.doi.org/10.1101/027201>.

Maier et al. 2015. Joint analysis of psychiatric disorders increases
accuracy of risk prediction for schizophrenia, bipolar disorder, and
major depressive disorder. Am J Hum Genet; 96(2):283-294.

Rodriguez-Alvarez, Maria Xose, et al. Correcting for spatial
heterogeneity in plant breeding experiments with P-splines. Spatial
Statistics 23 (2018): 52-71.

Searle. 1993. Applying the EM algorithm to calculating ML and REML
estimates of variance components. Paper invited for the 1993 American
Statistical Association Meeting, San Francisco.

Yu et al. 2006. A unified mixed-model method for association mapping
that accounts for multiple levels of relatedness. Genetics 38:203-208.

Tunnicliffe W. 1989. On the use of marginal likelihood in time series
model estimation. JRSS 51(1):15-27.
