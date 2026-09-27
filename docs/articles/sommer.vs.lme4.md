# Translating lme4 models to sommer

The sommer package was developed to provide R users with a powerful and
reliable multivariate mixed model solver for different genetic and
non-genetic analyses in diploid and polyploid organisms. This package
allows the user to estimate variance components for a mixed model with
the advantages of specifying the variance-covariance structure of the
random effects, specifying heterogeneous variances, and obtaining other
parameters such as BLUPs, BLUEs, residuals, fitted values, variances for
fixed and random effects, etc. The core algorithms of the package are
coded in C++ using the Armadillo library to optimize dense matrix
operations common in the derect-inversion algorithms. Although the
vignette shows examples using the mmes function with the default direct
inversion algorithm (henderson=FALSE) the Henderson’s approach can be
faster when the number of records surpasses the number of coefficients
to estimate and setting the henderson argument to TRUE can bring
significant speed ups.

The purpose of this vignette is to show how to translate the syntax
formula from `lme4` models to `sommer` models. Feel free to remove the
comment marks from the lme4 code so you can compare the results.

1.  Random slopes with same intercept
2.  Random slopes and random intercepts (without correlation)
3.  Random slopes and random intercepts (with correlation)
4.  Random slopes with a different intercept
5.  Other models not available in lme4

## 1) Random slopes

This is the simplest model people use when a random effect is desired
and the levels of the random effect are considered to have the same
intercept.

``` r
# install.packages("lme4")
# library(lme4)
library(sommer)
```

    ## Loading required package: Matrix

    ## Loading required package: MASS

    ## Loading required package: crayon

    ## Loading required package: enhancer

``` r
data(DT_sleepstudy, package="enhancer")
DT <- DT_sleepstudy
###########
## lme4
###########
# fm1 <- lmer(Reaction ~ Days + (1 | Subject), data=DT)
# summary(fm1) # or vc <- VarCorr(fm1); print(vc,comp=c("Variance"))
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  Subject  (Intercept) 1378.2   37.12   
#  Residual              960.5   30.99   
# Number of obs: 180, groups:  Subject, 18
###########
## sommer
###########
fm2 <- mmes(Reaction ~ Days,
            random= ~ Subject, 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
```

    ## Solver selected: ldlt

``` r
summary(fm2)$varcomp
```

    ##                term factor parameter  estimate StdError    Zratio
    ## 1 vsm(ism(Subject)) sigma2    sigma2 1377.9644 357.4925  3.854526
    ## 2   vsm(ism(units)) sigma2    sigma2  960.4978  75.6861 12.690544

## 2) Random slopes and random intercepts (without correlation)

This is the a model where you assume that the random effect has
different intercepts based on the levels of another variable. In
addition the `||` in `lme4` assumes that slopes and intercepts have no
correlation.

``` r
###########
## lme4
###########
# fm1 <- lmer(Reaction ~ Days + (Days || Subject), data=DT)
# summary(fm1) # or vc <- VarCorr(fm1); print(vc,comp=c("Variance"))
# Random effects:
#  Groups    Name        Variance Std.Dev.
#  Subject   (Intercept) 627.57   25.051  
#  Subject.1 Days         35.86    5.988  
#  Residual              653.58   25.565  
# Number of obs: 180, groups:  Subject, 18
###########
## sommer
###########
fm2 <- mmes(Reaction ~ Days,
            random= ~ Subject + vsm(ism(Days), ism(Subject)), 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
```

    ## Solver selected: ldlt

``` r
summary(fm2)$varcomp
```

    ##                           term   factor parameter  estimate  StdError    Zratio
    ## 1            vsm(ism(Subject))   sigma2    sigma2 627.55204 200.33887  3.132453
    ## 2 vsm(ism(Days), ism(Subject)) identity    sigma2  35.85471  10.26685  3.492281
    ## 3              vsm(ism(units))   sigma2    sigma2 653.58238  54.14683 12.070557

Notice that Days is a numerical (not factor) variable.

## 3) Random slopes and random intercepts (with correlation)

This is the a model where you assume that the random effect has
different intercepts based on the levels of another variable. In
addition a single `|` in `lme4` assumes that slopes and intercepts have
a correlation to be estimated.

``` r
###########
## lme4
###########
# fm1 <- lmer(Reaction ~ Days + (Days | Subject), data=DT)
# summary(fm1) # or # vc <- VarCorr(fm1); print(vc,comp=c("Variance"))
# Random effects:
#  Groups   Name        Variance Std.Dev. Corr
#  Subject  (Intercept) 612.10   24.741       
#           Days         35.07    5.922   0.07
#  Residual             654.94   25.592       
# Number of obs: 180, groups:  Subject, 18
###########
## sommer
###########
fm2 <- mmes(Reaction ~ Days, # henderson=TRUE,
            random= ~ vsm(ism(Days), ism(Subject)) , 
            nIters = 200, data=DT, tolParInv = 1e-6, verbose = FALSE)
```

    ## Solver selected: ldlt

``` r
summary(fm2)$varcomp
```

    ##                           term   factor parameter  estimate StdError    Zratio
    ## 1 vsm(ism(Days), ism(Subject)) identity    sigma2  52.70898 13.50174  3.903865
    ## 2              vsm(ism(units))   sigma2    sigma2 842.09199 66.36449 12.688894

``` r
cov2cor(fm2$theta[[1]])
```

    ##      [,1]
    ## [1,]    1

Notice that this last model require a new function called covm() which
creates the two random effects as before but now they have to be
encapsulated in covm() instead of just added.

## 4) Random slopes with a different intercept

This is the a model where you assume that the random effect has
different intercepts based on the levels of another variable but there’s
not a main effect. The 0 in the intercept in lme4 assumes that random
slopes interact with an intercept but without a main effect.

``` r
###########
## lme4
###########
# fm1 <- lmer(Reaction ~ Days + (0 + Days | Subject), data=DT)
# summary(fm1) # or vc <- VarCorr(fm1); print(vc,comp=c("Variance"))
# Random effects:
#  Groups   Name Variance Std.Dev.
#  Subject  Days  52.71    7.26   
#  Residual      842.03   29.02   
# Number of obs: 180, groups:  Subject, 18
###########
## sommer
###########
fm2 <- mmes(Reaction ~ Days,
            random= ~ vsm(ism(Days), ism(Subject)), 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
```

    ## Solver selected: ldlt

``` r
summary(fm2)$varcomp
```

    ##                           term   factor parameter  estimate StdError    Zratio
    ## 1 vsm(ism(Days), ism(Subject)) identity    sigma2  52.70898 13.50174  3.903865
    ## 2              vsm(ism(units))   sigma2    sigma2 842.09199 66.36449 12.688894

## 4) Other models available in sommer but not in lme4

One of the strengths of sommer is the availability of other variance
covariance structures. In this section we show 4 models available in
sommer that are not available in lme4 and might be useful.

``` r
library(orthopolynom)
## diagonal model
fm2 <- mmes(Reaction ~ Days,
            random= ~ vsm(dsm(Daysf), ism(Subject)), 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
```

    ## Solver selected: ldlt

``` r
summary(fm2)$varcomp
```

    ##                             term factor   parameter  estimate  StdError
    ## 1  vsm(dsm(Daysf), ism(Subject))   diag variance[0]  755.9399  296.6961
    ## 2  vsm(dsm(Daysf), ism(Subject))   diag variance[1]  812.7526  433.4803
    ## 3  vsm(dsm(Daysf), ism(Subject))   diag variance[2]  615.8352  407.7231
    ## 4  vsm(dsm(Daysf), ism(Subject))   diag variance[3] 1171.9460  490.2684
    ## 5  vsm(dsm(Daysf), ism(Subject))   diag variance[4] 1471.0567  543.7051
    ## 6  vsm(dsm(Daysf), ism(Subject))   diag variance[5] 2315.2861  710.5897
    ## 7  vsm(dsm(Daysf), ism(Subject))   diag variance[6] 3526.7361  971.6563
    ## 8  vsm(dsm(Daysf), ism(Subject))   diag variance[7] 2155.4792  679.2310
    ## 9  vsm(dsm(Daysf), ism(Subject))   diag variance[8] 3213.3728  905.5495
    ## 10 vsm(dsm(Daysf), ism(Subject))   diag variance[9] 4088.5269 1102.8136
    ## 11               vsm(ism(units)) sigma2      sigma2  263.7786  353.4608
    ##      Zratio
    ## 1  2.547860
    ## 2  1.874947
    ## 3  1.510425
    ## 4  2.390417
    ## 5  2.705615
    ## 6  3.258260
    ## 7  3.629613
    ## 8  3.173411
    ## 9  3.548534
    ## 10 3.707360
    ## 11 0.746274

``` r
## unstructured model
fm2 <- mmes(Reaction ~ Days,
            random= ~ vsm(usm(Daysf), ism(Subject)), 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
```

    ## Solver selected: ldlt

``` r
summary(fm2)$varcomp
```

    ##                             term factor       parameter   estimate   StdError
    ## 1  vsm(usm(Daysf), ism(Subject))     us     variance[0] 1259.35273 212.398210
    ## 2  vsm(usm(Daysf), ism(Subject))     us covariance[1,0] 1019.66534 203.692167
    ## 3  vsm(usm(Daysf), ism(Subject))     us covariance[2,0]  538.17404 149.364596
    ## 4  vsm(usm(Daysf), ism(Subject))     us covariance[3,0]  794.69035 188.512694
    ## 5  vsm(usm(Daysf), ism(Subject))     us covariance[4,0]  753.65790 191.189616
    ## 6  vsm(usm(Daysf), ism(Subject))     us covariance[5,0]  927.47612 245.954379
    ## 7  vsm(usm(Daysf), ism(Subject))     us covariance[6,0]  661.13179 291.872842
    ## 8  vsm(usm(Daysf), ism(Subject))     us covariance[7,0]  927.83239 260.538357
    ## 9  vsm(usm(Daysf), ism(Subject))     us covariance[8,0]  914.40531 275.188355
    ## 10 vsm(usm(Daysf), ism(Subject))     us covariance[9,0] 1433.53342 346.208352
    ## 11 vsm(usm(Daysf), ism(Subject))     us     variance[1] 1274.59032 257.742657
    ## 12 vsm(usm(Daysf), ism(Subject))     us covariance[2,1]  827.11046 184.755093
    ## 13 vsm(usm(Daysf), ism(Subject))     us covariance[3,1] 1135.06183 235.116929
    ## 14 vsm(usm(Daysf), ism(Subject))     us covariance[4,1] 1037.07500 226.919158
    ## 15 vsm(usm(Daysf), ism(Subject))     us covariance[5,1] 1176.35603 278.753024
    ## 16 vsm(usm(Daysf), ism(Subject))     us covariance[6,1]  848.17591 305.281008
    ## 17 vsm(usm(Daysf), ism(Subject))     us covariance[7,1]  923.79510 266.248407
    ## 18 vsm(usm(Daysf), ism(Subject))     us covariance[8,1] 1041.18556 291.412593
    ## 19 vsm(usm(Daysf), ism(Subject))     us covariance[9,1] 1509.70893 369.235916
    ## 20 vsm(usm(Daysf), ism(Subject))     us     variance[2]  858.82966 180.005389
    ## 21 vsm(usm(Daysf), ism(Subject))     us covariance[3,2] 1056.34756 214.542386
    ## 22 vsm(usm(Daysf), ism(Subject))     us covariance[4,2]  910.37612 198.322511
    ## 23 vsm(usm(Daysf), ism(Subject))     us covariance[5,2]  858.26734 224.666605
    ## 24 vsm(usm(Daysf), ism(Subject))     us covariance[6,2]  915.83771 273.154745
    ## 25 vsm(usm(Daysf), ism(Subject))     us covariance[7,2]  923.27624 241.640572
    ## 26 vsm(usm(Daysf), ism(Subject))     us covariance[8,2]  830.46628 243.779909
    ## 27 vsm(usm(Daysf), ism(Subject))     us covariance[9,2]  965.63016 288.081813
    ## 28 vsm(usm(Daysf), ism(Subject))     us     variance[3] 1618.37630 306.725744
    ## 29 vsm(usm(Daysf), ism(Subject))     us covariance[4,3] 1589.49115 304.097191
    ## 30 vsm(usm(Daysf), ism(Subject))     us covariance[5,3] 1671.36740 349.848796
    ## 31 vsm(usm(Daysf), ism(Subject))     us covariance[6,3] 1785.11239 419.663290
    ## 32 vsm(usm(Daysf), ism(Subject))     us covariance[7,3] 1280.44703 317.953384
    ## 33 vsm(usm(Daysf), ism(Subject))     us covariance[8,3] 1605.08053 369.793833
    ## 34 vsm(usm(Daysf), ism(Subject))     us covariance[9,3] 1738.71808 415.900452
    ## 35 vsm(usm(Daysf), ism(Subject))     us     variance[4] 1815.52612 346.362127
    ## 36 vsm(usm(Daysf), ism(Subject))     us covariance[5,4] 2001.65922 398.364102
    ## 37 vsm(usm(Daysf), ism(Subject))     us covariance[6,4] 2077.19657 465.258909
    ## 38 vsm(usm(Daysf), ism(Subject))     us covariance[7,4] 1547.79865 357.889197
    ## 39 vsm(usm(Daysf), ism(Subject))     us covariance[8,4] 2027.71310 427.501178
    ## 40 vsm(usm(Daysf), ism(Subject))     us covariance[9,4] 2194.80371 479.214379
    ## 41 vsm(usm(Daysf), ism(Subject))     us     variance[5] 2925.34776 569.932469
    ## 42 vsm(usm(Daysf), ism(Subject))     us covariance[6,5] 2602.55124 594.964314
    ## 43 vsm(usm(Daysf), ism(Subject))     us covariance[7,5] 1938.72246 461.233293
    ## 44 vsm(usm(Daysf), ism(Subject))     us covariance[8,5] 3056.09162 615.781640
    ## 45 vsm(usm(Daysf), ism(Subject))     us covariance[9,5] 3232.05273 672.271818
    ## 46 vsm(usm(Daysf), ism(Subject))     us     variance[6] 3983.31578 840.090155
    ## 47 vsm(usm(Daysf), ism(Subject))     us covariance[7,6] 2304.36475 563.670578
    ## 48 vsm(usm(Daysf), ism(Subject))     us covariance[8,6] 2926.95111 678.416579
    ## 49 vsm(usm(Daysf), ism(Subject))     us covariance[9,6] 2207.26114 663.712223
    ## 50 vsm(usm(Daysf), ism(Subject))     us     variance[7] 2526.21297 532.428706
    ## 51 vsm(usm(Daysf), ism(Subject))     us covariance[8,7] 2430.45766 551.210995
    ## 52 vsm(usm(Daysf), ism(Subject))     us covariance[9,7] 2394.44305 588.482206
    ## 53 vsm(usm(Daysf), ism(Subject))     us     variance[8] 3804.46447 758.495359
    ## 54 vsm(usm(Daysf), ism(Subject))     us covariance[9,8] 3842.02102 798.944176
    ## 55 vsm(usm(Daysf), ism(Subject))     us     variance[9] 4797.30728 968.474772
    ## 56               vsm(ism(units)) sigma2          sigma2   23.90661   7.800193
    ##      Zratio
    ## 1  5.929206
    ## 2  5.005913
    ## 3  3.603090
    ## 4  4.215580
    ## 5  3.941939
    ## 6  3.770927
    ## 7  2.265136
    ## 8  3.561212
    ## 9  3.322834
    ## 10 4.140667
    ## 11 4.945205
    ## 12 4.476794
    ## 13 4.827648
    ## 14 4.570240
    ## 15 4.220066
    ## 16 2.778345
    ## 17 3.469674
    ## 18 3.572891
    ## 19 4.088738
    ## 20 4.771133
    ## 21 4.923724
    ## 22 4.590382
    ## 23 3.820182
    ## 24 3.352816
    ## 25 3.820866
    ## 26 3.406623
    ## 27 3.351930
    ## 28 5.276298
    ## 29 5.226918
    ## 30 4.777399
    ## 31 4.253678
    ## 32 4.027153
    ## 33 4.340474
    ## 34 4.180611
    ## 35 5.241699
    ## 36 5.024698
    ## 37 4.464604
    ## 38 4.324798
    ## 39 4.743175
    ## 40 4.580004
    ## 41 5.132797
    ## 42 4.374298
    ## 43 4.203345
    ## 44 4.962947
    ## 45 4.807658
    ## 46 4.741534
    ## 47 4.088141
    ## 48 4.314386
    ## 49 3.325630
    ## 50 4.744697
    ## 51 4.409305
    ## 52 4.068845
    ## 53 5.015805
    ## 54 4.808873
    ## 55 4.953466
    ## 56 3.064874

``` r
## random regression (legendre polynomials)
fm2 <- mmes(Reaction ~ Days,
            random= ~ vsm(dsm(leg(Days,1)), ism(Subject)), 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
```

    ## Solver selected: ldlt

``` r
summary(fm2)$varcomp
```

    ##                                   term factor      parameter  estimate StdError
    ## 1 vsm(dsm(leg(Days, 1)), ism(Subject))   diag variance[leg0] 2815.1665 659.5459
    ## 2 vsm(dsm(leg(Days, 1)), ism(Subject))   diag variance[leg1]  473.6758 141.8338
    ## 3                      vsm(ism(units)) sigma2         sigma2  654.9374  54.5542
    ##      Zratio
    ## 1  4.268340
    ## 2  3.339654
    ## 3 12.005260

``` r
## unstructured random regression (legendre)
fm2 <- mmes(Reaction ~ Days,
            random= ~ vsm(usm(leg(Days,1)), ism(Subject)), 
            data=DT, tolParInv = 1e-6, verbose = FALSE)
```

    ## Solver selected: ldlt

``` r
summary(fm2)$varcomp
```

    ##                                   term factor             parameter  estimate
    ## 1 vsm(usm(leg(Days, 1)), ism(Subject))     us        variance[leg0] 2815.0498
    ## 2 vsm(usm(leg(Days, 1)), ism(Subject))     us covariance[leg1,leg0]  869.2572
    ## 3 vsm(usm(leg(Days, 1)), ism(Subject))     us        variance[leg1]  473.3989
    ## 4                      vsm(ism(units)) sigma2                sigma2  654.9239
    ##    StdError    Zratio
    ## 1 672.90016  4.183458
    ## 2 260.20880  3.340614
    ## 3 143.50367  3.298863
    ## 4  54.45607 12.026646

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
