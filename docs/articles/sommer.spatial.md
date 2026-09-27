# Spatial modeling using the sommer package

The sommer package was developed to provide R users with a powerful and
reliable multivariate mixed model solver for different genetic (in
diploid and polyploid organisms) and non-genetic analyses. This package
allows the user to estimate variance components in a mixed model with
the advantages of specifying the variance-covariance structure of the
random effects, specifying heterogeneous variances, and obtaining other
parameters such as BLUPs, BLUEs, residuals, fitted values, variances for
fixed and random effects, etc. The core algorithms of the package are
coded in C++ using the Armadillo library to optimize dense matrix
operations common in the derect-inversion algorithms.

This vignette is focused on showing the capabilities of sommer to fit
spatial models using the two dimensional splines models.

**SECTION 1: Introduction**

1.  Background in tensor products

**SECTION 2: Spatial models**

1.  Two dimensional splines (multiple spatial components)
2.  Two dimensional splines (single spatial component)
3.  Spatial models in multiple trials at once

## SECTION 1: Introduction

### Backgrounds in tensor products

TBD

## SECTION 2: Spatial models

### 1) Two dimensional splines (multiple spatial components)

In this example we show how to obtain the same results than using the
SpATS package. This is achieved by using the `spl2Db` function which is
a wrapper of the `tpsmmb` function.

``` r
library(sommer)
```

    ## Loading required package: Matrix

    ## Loading required package: MASS

    ## Loading required package: crayon

    ## Loading required package: enhancer

``` r
data(DT_yatesoats, package="enhancer")
DT <- DT_yatesoats
DT$row <- as.numeric(as.character(DT$row))
DT$col <- as.numeric(as.character(DT$col))
DT$R <- as.factor(DT$row)
DT$C <- as.factor(DT$col)

# SPATS MODEL
# m1.SpATS <- SpATS(response = "Y",
#                   spatial = ~ PSANOVA(col, row, nseg = c(14,21), degree = 3, pord = 2),
#                   genotype = "V", fixed = ~ 1,
#                   random = ~ R + C, data = DT,
#                   control = list(tolerance = 1e-04))
# 
# summary(m1.SpATS, which = "variances")
# 
# Spatial analysis of trials with splines 
# 
# Response:                   Y         
# Genotypes (as fixed):       V         
# Spatial:                    ~PSANOVA(col, row, nseg = c(14, 21), degree = 3, pord = 2)
# Fixed:                      ~1        
# Random:                     ~R + C    
# 
# 
# Number of observations:        72
# Number of missing data:        0
# Effective dimension:           17.09
# Deviance:                      483.405
# 
# Variance components:
#                   Variance            SD     log10(lambda)
# R                 1.277e+02     1.130e+01           0.49450
# C                 2.673e-05     5.170e-03           7.17366
# f(col)            4.018e-15     6.339e-08          16.99668
# f(row)            2.291e-10     1.514e-05          12.24059
# f(col):row        1.025e-04     1.012e-02           6.59013
# col:f(row)        8.789e+01     9.375e+00           0.65674
# f(col):f(row)     8.036e-04     2.835e-02           5.69565
# 
# Residual          3.987e+02     1.997e+01 

# SOMMER MODEL
M <- spl2Dmats(x.coord.name = "col", y.coord.name = "row", data=DT, 
               nseg =c(14,21), degree = c(3,3), penaltyord = c(2,2) 
               )
mix <- mmes(Y~V, henderson = TRUE,
            random=~ R + C + vsm(ism(M$fC)) + vsm(ism(M$fR)) + 
              vsm(ism(M$fC.R)) + vsm(ism(M$C.fR)) +
              vsm(ism(M$fC.fR)),
            rcov=~units, verbose=FALSE,
            data=M$data)
```

    ## Solver selected: cholmod

``` r
summary(mix)$varcomp
```

    ##                term factor parameter     estimate     StdError      Zratio
    ## 1       vsm(ism(R)) sigma2    sigma2 100.87679372  59.65910862 1.690886707
    ## 2       vsm(ism(C)) sigma2    sigma2 178.49151202 120.05530581 1.486744054
    ## 3    vsm(ism(M$fC)) sigma2    sigma2   0.57442858   3.31644824 0.173205954
    ## 4    vsm(ism(M$fR)) sigma2    sigma2   0.01765388   0.10192441 0.173205656
    ## 5  vsm(ism(M$fC.R)) sigma2    sigma2   0.01632156   0.09423189 0.173206290
    ## 6  vsm(ism(M$C.fR)) sigma2    sigma2   0.01585849  13.26911559 0.001195143
    ## 7 vsm(ism(M$fC.fR)) sigma2    sigma2   0.01484255   0.08569330 0.173205439
    ## 8   vsm(ism(units)) sigma2    sigma2 501.02215833  73.07425037 6.856343456

### 2) Two dimensional splines in single field (single spatial component)

To reduce the computational burden of fitting multiple spatial kernels
`sommer` provides a single spatial kernel method through the `spl2Da`
function. This as will be shown, can produce similar results to the more
flexible model. Use the one that fits better your needs.

``` r
# SOMMER MODEL
mix <- mmes(Y~V,
            random=~ R + C +
              vsm(ism(spl2Dc(row,col)$Z$`A:all`)),
            rcov=~units, verbose=FALSE,
            data=DT)
```

    ## Solver selected: cholmod

``` r
summary(mix)$varcomp
```

    ##                                   term factor parameter estimate  StdError
    ## 1                          vsm(ism(R)) sigma2    sigma2 112.0185  58.43723
    ## 2                          vsm(ism(C)) sigma2    sigma2 157.4612 114.33871
    ## 3 vsm(ism(spl2Dc(row, col)$Z$`A:all`)) sigma2    sigma2 406.7127 274.36989
    ## 4                      vsm(ism(units)) sigma2    sigma2 405.0701  67.19725
    ##     Zratio
    ## 1 1.916902
    ## 2 1.377147
    ## 3 1.482352
    ## 4 6.028075

### 3) Spatial models in multiple trials at once

Sometimes we want to fit heterogeneous variance components when e.g.,
have multiple trials or different locations. The spatial models can also
be fitted that way using the `at.var` and `at.levels` arguments. The
first argument expects a variable that will define the levels at which
the variance components will be fitted. The second argument is a way for
the user to specify the levels at which the spatial kernels should be
fitted if the user doesn’t want to fit it for all levels (e.g., trials
or fields).

``` r
DT2 <- rbind(DT,DT)
DT2$Y <- DT2$Y + rnorm(length(DT2$Y))
DT2$trial <- c(rep("A",nrow(DT)),rep("B",nrow(DT)))
head(DT2)
```

    ##   row col         Y   N          V  B         MP R C trial
    ## 1   1   1  89.59996 0.2    Victory B2    Victory 1 1     A
    ## 2   2   1  61.25532   0    Victory B2    Victory 2 1     A
    ## 3   3   1 118.56274 0.4 Marvellous B2 Marvellous 3 1     A
    ## 4   4   1 143.99443 0.6 Marvellous B2 Marvellous 4 1     A
    ## 5   5   1 149.62155 0.6 GoldenRain B2 GoldenRain 5 1     A
    ## 6   6   1 109.14841 0.2 GoldenRain B2 GoldenRain 6 1     A

``` r
# SOMMER MODEL
mix <- mmes(Y~V,
            random=~ R + C +
              vsm(dsm(trial),ism(spl2Dc(row,col)$Z$`A:all`)),
            rcov=~units, verbose=FALSE,
            data=DT2)
```

    ## Solver selected: cholmod

``` r
summary(mix)$varcomp
```

    ##                                               term factor   parameter estimate
    ## 1                                      vsm(ism(R)) sigma2      sigma2 191.6885
    ## 2                                      vsm(ism(C)) sigma2      sigma2 186.1152
    ## 3 vsm(dsm(trial), ism(spl2Dc(row, col)$Z$`A:all`))   diag variance[A] 262.1648
    ## 4 vsm(dsm(trial), ism(spl2Dc(row, col)$Z$`A:all`))   diag variance[B] 216.0412
    ## 5                                  vsm(ism(units)) sigma2      sigma2 342.4939
    ##    StdError   Zratio
    ## 1  59.62629 3.214831
    ## 2 115.96914 1.604868
    ## 3 164.58851 1.592850
    ## 4 190.33843 1.135037
    ## 5  37.62855 9.101970

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
