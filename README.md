# sommer: Solving Mixed Model Equations in R

Structural multivariate-univariate linear mixed model solver for estimation of multiple random effects with unknown variance-covariance structures (e.g., heterogeneous and unstructured) and known covariance among levels of random effects (e.g., pedigree and genomic relationship matrices) ([Covarrubias-Pazaran, 2016](https://doi.org/10.1371/journal.pone.0156744); [Maier et al., 2015](https://doi.org/10.1016/j.ajhg.2014.12.006); Jensen et al., 1997). REML estimates can be obtained using the Direct-Inversion Newton-Raphson and Direct-Inversion Average Information algorithms for the problems r x r (r being the number of records) or using the Henderson-based average information algorithm for the problem c x c (c being the number of coefficients to estimate). Spatial models can also be fitted using the two-dimensional spline functionality available.

## Installation

You can install the development version of `sommer` from GitHub:

``` r
remotes::install_github('covaruber/sommer')
```

## Using `sommer` code and the GPL license

Please be aware that the `sommer` project is open-source software released under the GNU General Public License (GPL). The GPL gives you the freedom to use, study, modify, and redistribute the software, subject to the conditions of the license.

In particular, if you copy, modify, or create a work based on GPL-covered `sommer` source code and **distribute that work**, the GPL imposes requirements on how the resulting work may be licensed and distributed. Among other things, recipients must retain the freedoms provided by the GPL, and when GPL-covered software is distributed in binary or other non-source form, the corresponding source code must also be made available in accordance with the GPL.

Therefore, the fact that `sommer` is free and open-source does **not** mean that its source code can be copied into another software package, modified, redistributed as a proprietary product, and kept secret. If your project contains or is based on GPL-covered `sommer` code, please carefully review the GPL requirements before distributing it.

The purpose of the `sommer` project is to encourage open collaboration, reproducibility, and continued improvement of statistical and quantitative-genetics software. Contributions and derivative work are welcome, but they should respect the freedoms and obligations established by the GPL.

Please read the full GPL license distributed with `sommer` before incorporating `sommer` source code into another project.

 
## Development

The sommer package is under active development. If you are an expert in mixed models, statistics or programming and you know how to implement of the following:

+ generalized linear models
+ ...

please help us to take sommer to the next level. Drop me an email or make a pull request through github :)  
