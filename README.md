# sommer: Solving Mixed Model Equations in R

Structural multivariate-univariate linear mixed model solver for estimation of multiple random effects with unknown variance-covariance structures (e.g., heterogeneous and unstructured) and known covariance among levels of random effects (e.g., pedigree and genomic relationship matrices) ([Covarrubias-Pazaran, 2016](https://doi.org/10.1371/journal.pone.0156744); [Maier et al., 2015](https://doi.org/10.1016/j.ajhg.2014.12.006); Jensen et al., 1997). REML estimates can be obtained using the Direct-Inversion Newton-Raphson and Direct-Inversion Average Information algorithms for the problems r x r (r being the number of records) or using the Henderson-based average information algorithm for the problem c x c (c being the number of coefficients to estimate). Spatial models can also be fitted using the two-dimensional spline functionality available.

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
