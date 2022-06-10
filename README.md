# sommer: Solving Mixed Model Equations in R

Structural multivariate-univariate linear mixed model solver for estimation of multiple random effects and unknown variance-covariance structures (i.e. heterogeneous and unstructured variance models) ([Covarrubias-Pazaran, 2016](https://doi.org/10.1371/journal.pone.0156744); [Maier et al., 2015](https://doi.org/10.1016/j.ajhg.2014.12.006)). REML estimates can be obtained using the Direct-Inversion Newton-Raphson and Direct-Inversion Average Information algorithms. Designed for genomic prediction and genome wide association studies (GWAS), particularly focused in the p > n problem (more coefficients than observations) and dense known covariance structures for levels of random effects. Spatial models can also be fitted using i.e. the two-dimensional spline functionality available in sommer.

## Installation

You can install the development version of `sommer` from GitHub:

``` r
devtools::install_github('covaruber/sommer')
```

## Vignettes

 - [Quick start for the sommer package](https://cran.r-project.org/package=sommer/vignettes/v1.sommer.quick.start.pdf)
 - [Moving to newer versions of sommer](https://cran.r-project.org/package=sommer/vignettes/v2.sommer.changes.and.faqs.pdf)
 - [Quantitative genetics using the sommer package](https://cran.r-project.org/package=sommer/vignettes/v3.sommer.qg.pdf)
 - [GxE models in sommer](https://cran.r-project.org/package=sommer/vignettes/v4.sommer.gxe.pdf)
 - [lme4 vs sommer](https://cran.r-project.org/package=sommer/vignettes/v5.sommer.vs.lme4.pdf)
 
## Development

The sommer package is under active development. If you are an expert in mixed models, statistics or programming and you know how to implement of the following:

+ the minimum degree ordering algorithm 
+ the symbolic cholesky factorization
+ factor analytic structure
+ generalized linear models

please help us to take sommer to the next level. Drop me an email or push some changes through github :)  
