# **m**ixed **m**odel **e**quations for **r** records

The mmer function uses the Direct-Inversion Newton-Raphson or Average
Information coded in C++ using the Armadillo library to optimize dense
matrix operations common in genomic selection models. These algorithms
are **intended to be used for problems of the type c \> r (more
coefficients to estimate than records in the dataset) and/or dense
matrices**. For problems with sparse data, or problems of the type r \>
c (more records in the dataset than coefficients to estimate), the
MME-based algorithm in the
[`mmes`](https://covaruber.github.io/sommer/reference/mmes.md) function
is faster and we recommend to shift to use that function.

## Usage

``` r
mmer(fixed, random, rcov, data, weights, W, nIters=20, tolParConvLL = 1e-03, 
     tolParInv = 1e-06, init=NULL, constraints=NULL,method="NR", getPEV=TRUE,
     naMethodX="exclude", naMethodY="exclude",returnParam=FALSE, 
     dateWarning=TRUE,date.warning=TRUE,verbose=TRUE, reshapeOutput=TRUE, stepWeight=NULL,
     emWeight=NULL, contrasts=NULL)
```

## Arguments

- fixed:

  A formula specifying the **response variable(s)** **and fixed
  effects**, e.g.:

  *response ~ covariate* for univariate models

  *cbind(response.i,response.j) ~ covariate* for multivariate models

  The `fcm` function can be used to constrain fixed effects in
  multi-response models.

- random:

  A formula specifying the name of the **random effects**, e.g. *random=
  ~ genotype + year*.

  Useful functions can be used to fit heterogeneous variances and other
  special models (*see 'Special Functions' in the Details section for
  more information*):

  [`vsr`](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md)`(...,Gu,Gti,Gtc)`
  is the main function to specify variance models and special structures
  for random effects. On the ... argument you provide the unknown
  variance-covariance structures (e.g., usr,dsr,atr,csr) and the random
  effect where such covariance structure will be used (the random effect
  of interest). Gu is used to provide known covariance matrices among
  the levels of the random effect, Gti initial values and Gtc for
  constraints. Auxiliar functions for building the variance models are:

  \*\*
  [`dsr`](https://covaruber.github.io/sommer/reference/dsr.md)`(x)`,
  [`usr`](https://covaruber.github.io/sommer/reference/usr.md)`(x)`,
  [`csr`](https://covaruber.github.io/sommer/reference/csr.md)`(x)` and
  [`atr`](https://covaruber.github.io/sommer/reference/atr.md)`(x,levs)`
  can be used to specify unknown diagonal, unstructured and customized
  unstructured and diagonal covariance structures to be estimated by
  REML.

  \*\*
  [`unsm`](https://covaruber.github.io/sommer/reference/unsm.md)`(x)`,
  [`fixm`](https://covaruber.github.io/sommer/reference/fixm.md)`(x)`
  and [`diag`](https://rdrr.io/r/base/diag.html)`(x)` can be used to
  build easily matrices to specify constraints in the Gtc argument of
  the
  [`vsr()`](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md)
  function.

  \*\* [`overlay()`](https://rdrr.io/pkg/enhancer/man/overlay.html),
  [`spl2Da()`](https://covaruber.github.io/sommer/reference/spl2Dc.md),
  [`spl2Db()`](https://covaruber.github.io/sommer/reference/spl2Dc.md),
  and [`leg()`](https://rdrr.io/pkg/enhancer/man/leg.html) functions can
  be used to specify overlayed of design matrices of random effects, two
  dimensional spline and random regression models within the
  [`vsr()`](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md)
  function.

  `gvsr(...,Gu,Guc,Gti,Gtc)` is an alternative function to specify
  general variance structures between different random effects. An
  special case in the indirect genetic effect models. Is similar to the
  vsr function but in the ... argument the different random effects are
  provided.

- rcov:

  A formula specifying the name of the **error term**, e.g., *rcov= ~
  units*.

  Special heterogeneous and special variance models and constraints for
  the residual part are the same used on the random term but the name of
  the random effect is always "units" which can be thought as a column
  with as many levels as rows in the data, e.g.,
  *rcov=~vsr(dsr(covariate),units)*

- data:

  A data frame containing the variables specified in the formulas for
  response, fixed, and random effects.

- weights:

  Name of the covariate for weights. To be used for the product R =
  Wsi\*R\*Wsi, where \* is the matrix product, Wsi is the square root of
  the inverse of W and R is the residual matrix.

- W:

  Alternatively, instead of providing a vector of weights the user can
  specify an entire W matrix (e.g., when covariances exist). To be used
  first to produce Wis = solve(chol(W)), and then calculate R =
  Wsi\*R\*Wsi.t(), where \* is the matrix product, and R is the residual
  matrix. Only one of the arguments weights or W should be used. If both
  are indicated W will be given the preference.

- nIters:

  Maximum number of iterations allowed.

- tolParConvLL:

  Convergence criteria for the change in log-likelihood.

- tolParInv:

  Tolerance parameter for matrix inverse used when singularities are
  encountered in the estimation procedure.

- init:

  Initial values for the variance components. By default this is NULL
  and initial values for the variance components are provided by the
  algorithm, but in case the user want to provide initial values for ALL
  var-cov components this argument is functional. It has to be provided
  as a list, where each list element corresponds to one random effect
  (1x1 matrix) and if multitrait model is pursued each element of the
  list is a matrix of variance covariance components among traits for
  such random effect. Initial values can also be provided in the Gti
  argument of the
  [vsr](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md)
  function. Is highly encouraged to use the Gti and Gtc arguments of the
  [vsr](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md)
  function instead of this argument, but these argument can be used to
  provide all initial values at once

- constraints:

  When initial values are provided these have to be accompanied by their
  constraints. See the
  [vsr](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md)
  function for more details on the constraints. Is highly encouraged to
  use the Gti and Gtc arguments of the
  [vsr](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md)
  function instead of this argument but these argument can be used to
  provide all constraints at once.

- method:

  This refers to the method or algorithm to be used for estimating
  variance components. Direct-inversion Newton-Raphson **NR** and
  Average Information **AI** (Tunnicliffe 1989; Gilmour et al. 1995; Lee
  et al. 2015).

- getPEV:

  A TRUE/FALSE value indicating if the program should return the
  predicted error variance and variance for random effects. This option
  is provided since this can take a long time for certain models where p
  is \> n by a big extent.

- naMethodX:

  One of the two possible values; "include" or "exclude". If "include"
  is selected then the function will impute the X matrices for fixed
  effects with the median value. If "exclude" is selected it will get
  rid of all rows with missing values for the X (fixed) covariates. The
  default is "exclude". The "include" option should be used carefully.

- naMethodY:

  One of the three possible values; "include", "include2" or "exclude"
  (default) to treat the observations in response variable to be used in
  the estimation of variance components. The first option "include" will
  impute the response variables for all rows with the median value,
  whereas "include2" imputes the responses only for rows where there is
  observation(s) for at least one of the responses (only available in
  the multi-response models). If "exclude" is selected (default) it will
  get rid of rows in response(s) where missing values are present for at
  least one of the responses.

- returnParam:

  A TRUE/FALSE value to indicate if the program should return the
  parameters to be used for fitting the model instead of fitting the
  model.

- dateWarning:

  A TRUE/FALSE value to indicate if the program should warn you when is
  time to update the sommer package.

- date.warning:

  A TRUE/FALSE value to indicate if the program should warn you when is
  time to update the sommer package. This argument will be removed soon,
  just left for backcompatibility.

- verbose:

  A TRUE/FALSE value to indicate if the program should return the
  progress of the iterative algorithm.

- reshapeOutput:

  A TRUE/FALSE value to indicate if the output should be reshaped to be
  easier to interpret for the user, some information is missing from the
  multivariate models for an easy interpretation.

- stepWeight:

  A vector of values (of length equal to the number of iterations)
  indicating the weight used to multiply the update (delta) for variance
  components at each iteration. If NULL the 1st iteration will be
  multiplied by 0.5, the 2nd by 0.7, and the rest by 0.9. This argument
  can help to avoid that variance components go outside the parameter
  space in the initial iterations which doesn't happen very often with
  the NR method but it can be detected by looking at the behavior of the
  likelihood. In that case you may want to give a smaller weight to the
  initial 8-10 iterations.

- emWeight:

  A vector of values (of length equal to the number of iterations)
  indicating with values between 0 and 1 the weight assigned to the EM
  information matrix. And the values 1 - emWeight will be applied to the
  NR/AI information matrix to produce a joint information matrix.

- contrasts:

  an optional list. See the contrasts.arg of model.matrix.default.

## Details

The use of this function requires a good understanding of mixed models.
Please review the 'sommer.quick.start' vignette and pay attention to
details like format of your random and fixed variables (e.g. character
and factor variables have different properties when returning BLUEs or
BLUPs, please see the 'sommer.changes.and.faqs' vignette).

**For tutorials** on how to perform different analysis with sommer
please look at the vignettes by typing in the terminal:

vignette("v1.sommer.quick.start")

vignette("v2.sommer.changes.and.faqs")

vignette("v3.sommer.qg")

vignette("v4.sommer.gxe")

**Citation**

Type *citation("sommer")* to know how to cite the sommer package in your
publications.

**Special variance structures**

[`vsr`](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md)`(`[`atr`](https://covaruber.github.io/sommer/reference/atr.md)`(x,levels),y)`

can be used to specify heterogeneous variance for the "y" covariate at
specific levels of the covariate "x", e.g.,
*random=~vsr(at(Location,c("A","B")),ID)* fits a variance component for
ID at levels A and B of the covariate Location.

[`vsr`](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md)`(`[`dsr`](https://covaruber.github.io/sommer/reference/dsr.md)`(x),y)`

can be used to specify a diagonal covariance structure for the "y"
covariate for all levels of the covariate "x", e.g.,
*random=~vsr(dsr(Location),ID)* fits a variance component for ID at all
levels of the covariate Location.

[`vsr`](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md)`(`[`usr`](https://covaruber.github.io/sommer/reference/usr.md)`(x),y)`

can be used to specify an unstructured covariance structure for the "y"
covariate for all levels of the covariate "x", e.g.,
*random=~vsr(usr(Location),ID)* fits variance and covariance components
for ID at all levels of the covariate Location.

[`vsr`](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md)`(`[`overlay`](https://rdrr.io/pkg/enhancer/man/overlay.html)`(...,rlist=NULL,prefix=NULL))`

can be used to specify overlay of design matrices between consecutive
random effects specified, e.g., *random=~vsr(overlay(male,female))*
overlays (overlaps) the incidence matrices for the male and female
random effects to obtain a single variance component for both effects.
The \`rlist\` argument is a list with each element being a numeric value
that multiplies the incidence matrix to be overlayed. See
[`overlay`](https://rdrr.io/pkg/enhancer/man/overlay.html) for
details.Can be combined with vsr().

[`vsr`](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md)`(`[`leg`](https://rdrr.io/pkg/enhancer/man/leg.html)`(x,n),y)`

can be used to fit a random regression model using a numerical variable
`x` that marks the trayectory for the random effect `y`. The leg
function can be combined with the special functions `dsr`, `usr` `at`
and `csr`. For example *random=~vsr(leg(x,1),y)* or
*random=~vsr(usr(leg(x,1)),y)*.

[`vsr`](https://covaruber.github.io/sommer/reference/DEPRECATED_VSR.md)`(x,Gtc=fcm(v))`

can be used to constrain fixed effects in the multi-response mixed
models. This is a vector that specifies if the fixed effect is to be
estimated for such trait. For example *fixed=cbind(response.i,
response.j)~vsr(Rowf, Gtc=fcm(c(1,0)))* means that the fixed effect Rowf
should only be estimated for the first response and the second should
only have the intercept.

`gvsr(x,y)`

can be used to fit variance and covariance parameters between two or
more random effects. For example, indirect genetic effect models.

[`spl2Da`](https://covaruber.github.io/sommer/reference/spl2Dc.md)`(x.coord, y.coord, at.var, at.levels))`

can be used to fit a 2-dimensional spline (e.g., spatial modeling) using
coordinates `x.coord` and `y.coord` (in numeric class) assuming a single
variance component. The 2D spline can be fitted at specific levels using
the `at.var` and `at.levels` arguments. For example
*random=~spl2Da(x.coord=Row,y.coord=Range,at.var=FIELD)*.

[`spl2Db`](https://covaruber.github.io/sommer/reference/spl2Dc.md)`(x.coord, y.coord, at.var, at.levels))`

can be used to fit a 2-dimensional spline (e.g., spatial modeling) using
coordinates `x.coord` and `y.coord` (in numeric class) assuming multiple
variance components. The 2D spline can be fitted at specific levels
using the `at.var` and `at.levels` arguments. For example
*random=~spl2Db(x.coord=Row,y.coord=Range,at.var=FIELD)*.

**S3 methods**

S3 methods are available for some parameter extraction such as
[`fitted.mmer`](https://covaruber.github.io/sommer/reference/fitted_mmes.md),
[`residuals.mmer`](https://covaruber.github.io/sommer/reference/residuals_mmes.md),
[`summary.mmer`](https://covaruber.github.io/sommer/reference/DEPRECATED_summary_mmer.md),
[`randef`](https://covaruber.github.io/sommer/reference/randef.md),
[`coef.mmer`](https://covaruber.github.io/sommer/reference/coef_mmes.md),
[`anova.mmer`](https://covaruber.github.io/sommer/reference/anova_mmes.md),
[`plot.mmer`](https://covaruber.github.io/sommer/reference/plot_mmes.md),
and
[`predict.mmer`](https://covaruber.github.io/sommer/reference/predict_mmes.md)
to obtain adjusted means. In addition, the
[`vpredict`](https://covaruber.github.io/sommer/reference/vpredict.md)
function (replacement of the pin function) can be used to estimate
standard errors for linear combinations of variance components (e.g.,
ratios like h2).

**Additional Functions**

Additional functions for genetic analysis have been included such as
relationship matrix building
([`A.mat`](https://covaruber.github.io/sommer/reference/A.mat.md),
[`D.mat`](https://covaruber.github.io/sommer/reference/D.mat.md),
[`E.mat`](https://covaruber.github.io/sommer/reference/E.mat.md),
[`H.mat`](https://covaruber.github.io/sommer/reference/H.mat.md)), build
a genotypic hybrid marker matrix
([`build.HMM`](https://rdrr.io/pkg/enhancer/man/build.HMM.html)), plot
of genetic maps
([`map.plot`](https://rdrr.io/pkg/enhancer/man/map.plot.html)), and
manhattan plots
([`manhattan`](https://rdrr.io/pkg/enhancer/man/manhattan.html)). If you
need to build a pedigree-based relationship matrix use the `getA`
function from the pedigreemm package.

**Bug report and contact**

If you have any technical questions or suggestions please post it in
https://stackoverflow.com or https://stats.stackexchange.com

If you have any bug report please go to
https://github.com/covaruber/sommer or send me an email to address it
asap, just make sure you have read the vignettes carefully before
sending your question.

**Example Datasets**

The package has been equiped with several datasets to learn how to use
the sommer package:

\*
[`DT_halfdiallel`](https://rdrr.io/pkg/enhancer/man/DT_halfdiallel.html),
[`DT_fulldiallel`](https://rdrr.io/pkg/enhancer/man/DT_fulldiallel.html)
and [`DT_mohring`](https://rdrr.io/pkg/enhancer/man/DT_mohring.html)
datasets have examples to fit half and full diallel designs.

\* [`DT_h2`](https://rdrr.io/pkg/enhancer/man/DT_h2.html) to calculate
heritability

\*
[`DT_cornhybrids`](https://rdrr.io/pkg/enhancer/man/DT_cornhybrids.html)
and [`DT_technow`](https://rdrr.io/pkg/enhancer/man/DT_technow.html)
datasets to perform genomic prediction in hybrid single crosses

\* [`DT_wheat`](https://rdrr.io/pkg/enhancer/man/DT_wheat.html) dataset
to do genomic prediction in single crosses in species displaying only
additive effects.

\* [`DT_cpdata`](https://rdrr.io/pkg/enhancer/man/DT_cpdata.html)
dataset to fit genomic prediction models within a biparental population
coming from 2 highly heterozygous parents including additive, dominance
and epistatic effects.

\* [`DT_polyploid`](https://rdrr.io/pkg/enhancer/man/DT_polyploid.html)
to fit genomic prediction and GWAS analysis in polyploids.

\* [`DT_gryphon`](https://rdrr.io/pkg/enhancer/man/DT_gryphon.html) data
contains an example of an animal model including pedigree information.

\* [`DT_btdata`](https://rdrr.io/pkg/enhancer/man/DT_btdata.html)
dataset contains an animal (birds) model.

\* [`DT_legendre`](https://rdrr.io/pkg/enhancer/man/DT_legendre.html)
simulated dataset for random regression model.

\*
[`DT_sleepstudy`](https://rdrr.io/pkg/enhancer/man/DT_sleepstudy.html)
dataset to know how to translate lme4 models to sommer models.

\* [`DT_ige`](https://rdrr.io/pkg/enhancer/man/DT_ige.html) dataset to
show how to fit indirect genetic effect models.

**Models Enabled**

For details about the models enabled and more information about the
covariance structures please check the help page of the package
([`sommer`](https://covaruber.github.io/sommer/reference/sommer-package.md)).

## Value

If all parameters are correctly indicated the program will return a list
with the following information:

- Vi:

  the inverse of the phenotypic variance matrix V^- = (ZGZ+R)^-1

- P:

  the projection matrix Vi - \[Vi\*(X\*Vi\*X)^-\*Vi\]

- sigma:

  a list with the values of the variance-covariance components with one
  list element for each random effect.

- sigma_scaled:

  a list with the values of the scaled variance-covariance components
  with one list element for each random effect.

- sigmaSE:

  Hessian matrix containing the variance-covariance for the variance
  components. SE's can be obtained taking the square root of the
  diagonal values of the Hessian.

- Beta:

  a data frame for trait BLUEs (fixed effects).

- VarBeta:

  a variance-covariance matrix for trait BLUEs

- U:

  a list (one element for each random effect) with a data frame for
  trait BLUPs.

- VarU:

  a list (one element for each random effect) with the
  variance-covariance matrix for trait BLUPs.

- PevU:

  a list (one element for each random effect) with the predicted error
  variance matrix for trait BLUPs.

- fitted:

  Fitted values y.hat=XB

- residuals:

  Residual values e = Y - XB

- AIC:

  Akaike information criterion

- BIC:

  Bayesian information criterion

- convergence:

  a TRUE/FALSE statement indicating if the model converged.

- monitor:

  The values of log-likelihood and variance-covariance components across
  iterations during the REML estimation.

- percChange:

  The percent change of variance components across iterations. There
  should be one column less than the number of iterations. Calculated as
  percChange = ((x_i/x_i-1) - 1) \* 100 where i is the ith iteration.

- dL:

  The vector of first derivatives of the likelihood with respect to the
  ith variance-covariance component.

- dL2:

  The matrix of second derivatives of the likelihood with respect to the
  i.j th variance-covariance component.

- method:

  The method for extimation of variance components specified by the
  user.

- call:

  Formula for fixed, random and rcov used.

- constraints:

  contraints used in the mixed models for the random effects.

- constraintsF:

  contraints used in the mixed models for the fixed effects.

- data:

  The dataset used in the model after removing missing records for the
  response variable.

- dataOriginal:

  The original dataset used in the model.

- terms:

  The name of terms for responses, fixed, random and residual effects in
  the model.

- termsN:

  The number of effects associated to fixed, random and residual effects
  in the model.

- sigmaVector:

  a vectorized version of the sigma element (variance-covariance
  components) to match easily the standard errors of the var-cov
  components stored in the element sigmaSE.

- reshapeOutput:

  The value provided to the mmer function for the argument with the same
  name.

## References

Covarrubias-Pazaran G. Genome assisted prediction of quantitative traits
using the R package sommer. PLoS ONE 2016, 11(6):
doi:10.1371/journal.pone.0156744

Covarrubias-Pazaran G. 2018. Software update: Moving the R package
sommer to multivariate mixed models for genome-assisted prediction. doi:
https://doi.org/10.1101/354639

Bernardo Rex. 2010. Breeding for quantitative traits in plants. Second
edition. Stemma Press. 390 pp.

Gilmour et al. 1995. Average Information REML: An efficient algorithm
for variance parameter estimation in linear mixed models. Biometrics
51(4):1440-1450.

Kang et al. 2008. Efficient control of population structure in model
organism association mapping. Genetics 178:1709-1723.

Lee, D.-J., Durban, M., and Eilers, P.H.C. (2013). Efficient
two-dimensional smoothing with P-spline ANOVA mixed models and nested
bases. Computational Statistics and Data Analysis, 61, 22 - 37.

Lee et al. 2015. MTG2: An efficient algorithm for multivariate linear
mixed model analysis based on genomic information. Cold Spring Harbor.
doi: http://dx.doi.org/10.1101/027201.

Maier et al. 2015. Joint analysis of psychiatric disorders increases
accuracy of risk prediction for schizophrenia, bipolar disorder, and
major depressive disorder. Am J Hum Genet; 96(2):283-294.

Rodriguez-Alvarez, Maria Xose, et al. Correcting for spatial
heterogeneity in plant breeding experiments with P-splines. Spatial
Statistics 23 (2018): 52-71.

Searle. 1993. Applying the EM algorithm to calculating ML and REML
estimates of variance components. Paper invited for the 1993 American
Statistical Association Meeting, San Francisco.

Yu et al. 2006. A unified mixed-model method for association mapping
that accounts for multiple levels of relatedness. Genetics 38:203-208.

Tunnicliffe W. 1989. On the use of marginal likelihood in time series
model estimation. JRSS 51(1):15-27.

Zhang et al. 2010. Mixed linear model approach adapted for genome-wide
association studies. Nat. Genet. 42:355-360.

## Author

Giovanny Covarrubias-Pazaran

## Examples

``` r
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
####=========================================####

####=========================================####
#### EXAMPLES
#### Different models with sommer
####=========================================####

data(DT_example, package="enhancer")
DT <- DT_example
head(DT)
#>                   Name     Env Loc Year     Block Yield    Weight
#> 33  Manistee(MSL292-A) CA.2013  CA 2013 CA.2013.1     4 -1.904711
#> 65          CO02024-9W CA.2013  CA 2013 CA.2013.1     5 -1.446958
#> 66  Manistee(MSL292-A) CA.2013  CA 2013 CA.2013.2     5 -1.516271
#> 67            MSL007-B CA.2011  CA 2011 CA.2011.2     5 -1.435510
#> 68           MSR169-8Y CA.2013  CA 2013 CA.2013.1     5 -1.469051
#> 103         AC05153-1W CA.2013  CA 2013 CA.2013.1     6 -1.307167

####=========================================####
#### Univariate homogeneous variance models  ####
####=========================================####

## Compound simmetry (CS) model
ans1 <- mmer(Yield~Env,
             random= ~ Name + Env:Name,
             rcov= ~ units,
             data=DT)
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -31.2668   21:17:43      0           0
#>     2      -23.2804   21:17:43      0           0
#>     3      -20.4746   21:17:43      0           0
#>     4      -20.1501   21:17:43      0           0
#>     5      -20.1454   21:17:43      0           0
#>     6      -20.1454   21:17:43      0           0
summary(ans1)
#> $groups
#>          Yield
#> Name        41
#> Env:Name   123
#> 
#> $varcomp
#>                       VarComp VarCompSE   Zratio Constraint
#> Name.Yield-Yield     3.681877 1.6909561 2.177394   Positive
#> Env:Name.Yield-Yield 5.173062 1.4952313 3.459707   Positive
#> units.Yield-Yield    4.366285 0.6470458 6.748031   Positive
#> 
#> $betas
#>   Trait      Effect  Estimate Std.Error   t.value
#> 1 Yield (Intercept) 16.496350 0.6854966 24.064816
#> 2 Yield  EnvCA.2012 -5.776758 0.7558134 -7.643101
#> 3 Yield  EnvCA.2013 -6.380478 0.7960468 -8.015204
#> 
#> $method
#> [1] "NR"
#> 
#> $logo
#>          logLik      AIC      BIC Method Converge
#> Value -20.14538 46.29075 55.95182     NR     TRUE
#> 
#> attr(,"class")
#> [1] "summary.mmer" "list"        

####===========================================####
#### Univariate heterogeneous variance models  ####
####===========================================####

## Compound simmetry (CS) + Diagonal (DIAG) model
ans2 <- mmer(Yield~Env,
             random= ~Name + vsr(dsr(Env),Name),
             rcov= ~ vsr(dsr(Env),units),
             data=DT)
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -31.2668   21:17:43      0           0
#>     2      -19.8549   21:17:43      0           0
#>     3      -15.9797   21:17:43      0           0
#>     4      -15.4374   21:17:43      0           0
#>     5      -15.43   21:17:43      0           0
#>     6      -15.4298   21:17:43      0           0
summary(ans2)
#> $groups
#>              Yield
#> Name            41
#> CA.2011:Name    41
#> CA.2012:Name    41
#> CA.2013:Name    41
#> 
#> $varcomp
#>                             VarComp VarCompSE   Zratio Constraint
#> Name.Yield-Yield           2.962851 1.4962000 1.980251   Positive
#> CA.2011:Name.Yield-Yield  10.146369 4.5073271 2.251083   Positive
#> CA.2012:Name.Yield-Yield   1.877530 1.8697568 1.004158   Positive
#> CA.2013:Name.Yield-Yield   6.629152 2.5028114 2.648682   Positive
#> CA.2011:units.Yield-Yield  4.942450 1.5245057 3.242001   Positive
#> CA.2012:units.Yield-Yield  5.724963 1.3123015 4.362536   Positive
#> CA.2013:units.Yield-Yield  2.559880 0.6399685 4.000010   Positive
#> 
#> $betas
#>   Trait      Effect  Estimate Std.Error   t.value
#> 1 Yield (Intercept) 16.507650 0.8268377 19.964800
#> 2 Yield  EnvCA.2012 -5.816861 0.8575288 -6.783284
#> 3 Yield  EnvCA.2013 -6.412375 0.9356112 -6.853674
#> 
#> $method
#> [1] "NR"
#> 
#> $logo
#>          logLik      AIC      BIC Method Converge
#> Value -15.42983 36.85965 46.52072     NR     TRUE
#> 
#> attr(,"class")
#> [1] "summary.mmer" "list"        

####===========================================####
####  Univariate unstructured variance models  ####
####===========================================####

ans3 <- mmer(Yield~Env,
             random=~ vsr(usr(Env),Name),
             rcov=~vsr(dsr(Env),units), 
             data=DT)
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -37.9059   21:17:43      0           0
#>     2      -17.9745   21:17:43      0           0
#>     3      -12.2427   21:17:43      0           0
#>     4      -11.5121   21:17:43      0           0
#>     5      -11.5001   21:17:43      0           0
#>     6      -11.4997   21:17:43      0           0
summary(ans3)
#> $groups
#>                      Yield
#> CA.2011:Name            41
#> CA.2012:CA.2011:Name    82
#> CA.2012:Name            41
#> CA.2013:CA.2011:Name    82
#> CA.2013:CA.2012:Name    82
#> CA.2013:Name            41
#> 
#> $varcomp
#>                                     VarComp VarCompSE    Zratio Constraint
#> CA.2011:Name.Yield-Yield         15.6650010 5.4206906 2.8898534   Positive
#> CA.2012:CA.2011:Name.Yield-Yield  6.1101600 2.4850649 2.4587527   Unconstr
#> CA.2012:Name.Yield-Yield          4.5296090 1.8208107 2.4876881   Positive
#> CA.2013:CA.2011:Name.Yield-Yield  6.3844808 3.0658977 2.0824181   Unconstr
#> CA.2013:CA.2012:Name.Yield-Yield  0.3929997 1.5233985 0.2579757   Unconstr
#> CA.2013:Name.Yield-Yield          8.5972750 2.4837814 3.4613654   Positive
#> CA.2011:units.Yield-Yield         4.9698460 1.5322540 3.2434870   Positive
#> CA.2012:units.Yield-Yield         5.6729333 1.3007862 4.3611574   Positive
#> CA.2013:units.Yield-Yield         2.5570940 0.6392821 3.9999462   Positive
#> 
#> $betas
#>   Trait      Effect  Estimate Std.Error   t.value
#> 1 Yield (Intercept) 16.331246 0.8137223 20.069802
#> 2 Yield  EnvCA.2012 -5.695837 0.7403976 -7.692944
#> 3 Yield  EnvCA.2013 -6.271081 0.8190942 -7.656117
#> 
#> $method
#> [1] "NR"
#> 
#> $logo
#>          logLik      AIC      BIC Method Converge
#> Value -11.49971 28.99943 38.66049     NR     TRUE
#> 
#> attr(,"class")
#> [1] "summary.mmer" "list"        

# \donttest{

####==========================================####
#### Multivariate homogeneous variance models ####
####==========================================####

## Multivariate Compound simmetry (CS) model
DT$EnvName <- paste(DT$Env,DT$Name)
ans4 <- mmer(cbind(Yield, Weight) ~ Env,
              random= ~ vsr(Name, Gtc = unsm(2)) + vsr(EnvName,Gtc = unsm(2)),
              rcov= ~ vsr(units, Gtc = unsm(2)),
              data=DT)
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      66.0395   21:17:43      0           0
#>     2      131.529   21:17:44      1           0
#>     3      162.769   21:17:44      1           0
#>     4      166.983   21:17:44      1           0
#>     5      167.025   21:17:44      1           0
#>     6      167.025   21:17:44      1           0
summary(ans4)
#> $groups
#>           Yield Weight
#> u:Name       41     41
#> u:EnvName    94     94
#> 
#> $varcomp
#>                           VarComp  VarCompSE   Zratio Constraint
#> u:Name.Yield-Yield      3.7089322 1.68116825 2.206164   Positive
#> u:Name.Yield-Weight     0.9070677 0.37944464 2.390514   Unconstr
#> u:Name.Weight-Weight    0.2243447 0.08774729 2.556714   Positive
#> u:EnvName.Yield-Yield   5.0920817 1.47878716 3.443418   Positive
#> u:EnvName.Yield-Weight  1.0268613 0.30766632 3.337581   Unconstr
#> u:EnvName.Weight-Weight 0.2100748 0.06661018 3.153795   Positive
#> u:units.Yield-Yield     4.3836757 0.64941305 6.750212   Positive
#> u:units.Yield-Weight    0.9077303 0.14144941 6.417349   Unconstr
#> u:units.Weight-Weight   0.2280088 0.03377178 6.751458   Positive
#> 
#> $betas
#>    Trait      Effect   Estimate Std.Error   t.value
#> 1  Yield (Intercept) 16.4093261 0.6783123 24.191401
#> 2 Weight (Intercept)  0.9805532 0.1497069  6.549820
#> 3  Yield  EnvCA.2012 -5.6844175 0.7474015 -7.605574
#> 4 Weight  EnvCA.2012 -1.1846280 0.1592542 -7.438600
#> 5  Yield  EnvCA.2013 -6.2951783 0.7850422 -8.018904
#> 6 Weight  EnvCA.2013 -1.3558524 0.1681056 -8.065483
#> 
#> $method
#> [1] "NR"
#> 
#> $logo
#>         logLik       AIC       BIC Method Converge
#> Value 167.0252 -322.0505 -298.5695     NR     TRUE
#> 
#> attr(,"class")
#> [1] "summary.mmer" "list"        

####=============================================####
#### Multivariate heterogeneous variance models  ####
####=============================================####

## Multivariate Compound simmetry (CS) + Diagonal (DIAG) model
ans5 <- mmer(cbind(Yield, Weight) ~ Env,
              random= ~ vsr(Name, Gtc = unsm(2)) + vsr(dsr(Env),Name, Gtc = unsm(2)),
              rcov= ~ vsr(dsr(Env),units, Gtc = unsm(2)),
              data=DT)
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      66.0395   21:17:45      1           0
#>     2      138.617   21:17:45      1           0
#>     3      172.682   21:17:45      1           0
#>     4      177.662   21:17:46      2           0
#>     5      177.801   21:17:46      2           0
#>     6      177.813   21:17:46      2           0
#>     7      177.815   21:17:47      3           0
#>     8      177.815   21:17:47      3           0
summary(ans5)
#> $groups
#>              Yield Weight
#> u:Name          41     41
#> CA.2011:Name    41     41
#> CA.2012:Name    41     41
#> CA.2013:Name    41     41
#> 
#> $varcomp
#>                                VarComp  VarCompSE    Zratio Constraint
#> u:Name.Yield-Yield          3.31936097 1.45268532 2.2849828   Positive
#> u:Name.Yield-Weight         0.79393338 0.32621238 2.4337929   Unconstr
#> u:Name.Weight-Weight        0.19085145 0.07502697 2.5437712   Positive
#> CA.2011:Name.Yield-Yield    8.70656990 4.01470472 2.1686701   Positive
#> CA.2011:Name.Yield-Weight   1.77892286 0.83926393 2.1196227   Unconstr
#> CA.2011:Name.Weight-Weight  0.35965596 0.17902862 2.0089300   Positive
#> CA.2012:Name.Yield-Yield    2.57109159 1.94951063 1.3188395   Positive
#> CA.2012:Name.Yield-Weight   0.33245253 0.39840053 0.8344681   Unconstr
#> CA.2012:Name.Weight-Weight  0.03841851 0.08595438 0.4469639   Positive
#> CA.2013:Name.Yield-Yield    5.46908017 2.16307290 2.5283846   Positive
#> CA.2013:Name.Yield-Weight   1.34713066 0.50478780 2.6687069   Unconstr
#> CA.2013:Name.Weight-Weight  0.32902387 0.12207588 2.6952406   Positive
#> CA.2011:units.Yield-Yield   4.93852081 1.52317838 3.2422472   Positive
#> CA.2011:units.Yield-Weight  0.99446886 0.32150226 3.0931940   Unconstr
#> CA.2011:units.Weight-Weight 0.23982305 0.07394418 3.2432986   Positive
#> CA.2012:units.Yield-Yield   5.73887436 1.31532905 4.3630712   Positive
#> CA.2012:units.Yield-Weight  1.28009125 0.30156568 4.2448174   Unconstr
#> CA.2012:units.Weight-Weight 0.31806306 0.07286389 4.3651675   Positive
#> CA.2013:units.Yield-Yield   2.56126578 0.63992925 4.0024202   Positive
#> CA.2013:units.Yield-Weight  0.44568736 0.12644887 3.5246449   Unconstr
#> CA.2013:units.Weight-Weight 0.12231798 0.03057261 4.0009007   Positive
#> 
#> $betas
#>    Trait      Effect   Estimate Std.Error   t.value
#> 1  Yield (Intercept) 16.4242562 0.7890504 20.815218
#> 2 Weight (Intercept)  0.9866188 0.1682863  5.862738
#> 3  Yield  EnvCA.2012 -5.7338723 0.8265664 -6.936977
#> 4 Weight  EnvCA.2012 -1.1998497 0.1698017 -7.066181
#> 5  Yield  EnvCA.2013 -6.3128448 0.8756733 -7.209132
#> 6 Weight  EnvCA.2013 -1.3620901 0.1914535 -7.114470
#> 
#> $method
#> [1] "NR"
#> 
#> $logo
#>         logLik       AIC       BIC Method Converge
#> Value 177.8154 -343.6308 -320.1497     NR     TRUE
#> 
#> attr(,"class")
#> [1] "summary.mmer" "list"        

####===========================================####
#### Multivariate unstructured variance models ####
####===========================================####

ans6 <- mmer(cbind(Yield, Weight) ~ Env,
              random= ~ vsr(usr(Env),Name, Gtc = unsm(2)),
              rcov= ~ vsr(dsr(Env),units, Gtc = unsm(2)),
              data=DT)
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      56.6189   21:17:48      1           0
#>     2      140.894   21:17:48      1           0
#>     3      176.238   21:17:49      2           0
#>     4      181.462   21:17:49      2           0
#>     5      181.688   21:17:50      3           0
#>     6      181.746   21:17:50      3           0
#>     7      181.77   21:17:50      3           0
#>     8      181.781   21:17:51      4           0
#>     9      181.787   21:17:51      4           0
#>     10      181.791   21:17:52      5           0
#>     11      181.793   21:17:52      5           0
#>     12      181.794   21:17:53      6           0
#>     13      181.794   21:17:53      6           0
summary(ans6)
#> $groups
#>                      Yield Weight
#> CA.2011:Name            41     41
#> CA.2012:CA.2011:Name    82     82
#> CA.2012:Name            41     41
#> CA.2013:CA.2011:Name    82     82
#> CA.2013:CA.2012:Name    82     82
#> CA.2013:Name            41     41
#> 
#> $varcomp
#>                                       VarComp  VarCompSE   Zratio Constraint
#> CA.2011:Name.Yield-Yield           15.6450590 5.35692033 2.920532   Positive
#> CA.2011:Name.Yield-Weight           3.3585731 1.14633483 2.929836   Unconstr
#> CA.2011:Name.Weight-Weight          0.7181736 0.24870736 2.887625   Positive
#> CA.2012:CA.2011:Name.Yield-Yield    6.5289181 2.48614965 2.626116   Positive
#> CA.2012:CA.2011:Name.Yield-Weight   1.3505115 0.52387542 2.577925   Unconstr
#> CA.2012:CA.2011:Name.Weight-Weight  0.2842006 0.11258722 2.524271   Positive
#> CA.2012:Name.Yield-Yield            4.7893181 1.86183450 2.572365   Positive
#> CA.2012:Name.Yield-Weight           0.8640033 0.38376793 2.251369   Unconstr
#> CA.2012:Name.Weight-Weight          0.1693101 0.08354311 2.026620   Positive
#> CA.2013:CA.2011:Name.Yield-Yield    5.9933574 2.93830062 2.039736   Positive
#> CA.2013:CA.2011:Name.Yield-Weight   1.4231956 0.64973332 2.190430   Unconstr
#> CA.2013:CA.2011:Name.Weight-Weight  0.3379104 0.14680112 2.301824   Positive
#> CA.2013:CA.2012:Name.Yield-Yield    2.0986982 1.44033972 1.457086   Positive
#> CA.2013:CA.2012:Name.Yield-Weight   0.5239908 0.32356125 1.619449   Unconstr
#> CA.2013:CA.2012:Name.Weight-Weight  0.1341957 0.07571863 1.772294   Positive
#> CA.2013:Name.Yield-Yield            8.6256730 2.47810836 3.480749   Positive
#> CA.2013:Name.Yield-Weight           2.1047727 0.58747655 3.582735   Unconstr
#> CA.2013:Name.Weight-Weight          0.5125436 0.14284861 3.588020   Positive
#> CA.2011:units.Yield-Yield           4.9515911 1.52694110 3.242817   Positive
#> CA.2011:units.Yield-Weight          0.9992924 0.32285607 3.095164   Unconstr
#> CA.2011:units.Weight-Weight         0.2410810 0.07431937 3.243851   Positive
#> CA.2012:units.Yield-Yield           5.7790208 1.32422902 4.364064   Positive
#> CA.2012:units.Yield-Weight          1.2913597 0.30407913 4.246789   Unconstr
#> CA.2012:units.Weight-Weight         0.3211785 0.07355847 4.366303   Positive
#> CA.2013:units.Yield-Yield           2.5566754 0.63882793 4.002135   Positive
#> CA.2013:units.Yield-Weight          0.4451622 0.12631251 3.524292   Unconstr
#> CA.2013:units.Weight-Weight         0.1222735 0.03056234 4.000791   Positive
#> 
#> $betas
#>    Trait      Effect   Estimate Std.Error   t.value
#> 1  Yield (Intercept) 16.3341888 0.8253829 19.789833
#> 2 Weight (Intercept)  0.9677418 0.1770436  5.466120
#> 3  Yield  EnvCA.2012 -5.6637309 0.7448686 -7.603664
#> 4 Weight  EnvCA.2012 -1.1855499 0.1604319 -7.389739
#> 5  Yield  EnvCA.2013 -6.2153229 0.8339869 -7.452543
#> 6 Weight  EnvCA.2013 -1.3406465 0.1805560 -7.425099
#> 
#> $method
#> [1] "NR"
#> 
#> $logo
#>         logLik       AIC       BIC Method Converge
#> Value 181.7945 -351.5889 -328.1079     NR     TRUE
#> 
#> attr(,"class")
#> [1] "summary.mmer" "list"        

####=========================================####
####=========================================####
#### EXAMPLE SET 2
#### 2 variance components
#### one random effect with variance covariance structure
####=========================================####
####=========================================####

data("DT_cpdata", package="enhancer")
DT <- DT_cpdata
GT <- GT_cpdata
MP <- MP_cpdata
head(DT)
#>        id Row Col Year      color  Yield FruitAver Firmness Rowf Colf
#> P003 P003   3   1 2014 0.10075269 154.67     41.93  588.917    3    1
#> P004 P004   4   1 2014 0.13891940 186.77     58.79  640.031    4    1
#> P005 P005   5   1 2014 0.08681502  80.21     48.16  671.523    5    1
#> P006 P006   6   1 2014 0.13408561 202.96     48.24  687.172    6    1
#> P007 P007   7   1 2014 0.13519278 174.74     45.83  601.322    7    1
#> P008 P008   8   1 2014 0.17406685 194.16     44.63  656.379    8    1
GT[1:4,1:4]
#>      scaffold_50439_2381 scaffold_39344_153 uneak_3436043 uneak_2632033
#> P003                   0                  0             0             1
#> P004                   0                  0             0             1
#> P005                   0                 -1             0             1
#> P006                  -1                 -1            -1             0
#### create the variance-covariance matrix
A <- A.mat(GT)
#### look at the data and fit the model
mix1 <- mmer(Yield~1,
             random=~vsr(id, Gu=A) + Rowf,
             rcov=~units,
             data=DT)
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -153.728   21:17:54      1           0
#>     2      -153.608   21:17:54      1           0
#>     3      -153.575   21:17:54      1           0
#>     4      -153.572   21:17:54      1           0
#>     5      -153.572   21:17:54      1           0
summary(mix1)$varcomp
#>                     VarComp VarCompSE    Zratio Constraint
#> u:id.Yield-Yield   708.9916  307.9387  2.302379   Positive
#> Rowf.Yield-Yield   839.7109  395.6582  2.122314   Positive
#> units.Yield-Yield 3238.5180  286.4714 11.304856   Positive

#### multi trait example
mix2 <- mmer(cbind(Yield,color)~1,
              random = ~ vsr(id, Gu=A, Gtc = unsm(2)) + # unstructured at trait level
                            vsr(Rowf, Gtc=diag(2)) + # diagonal structure at trait level
                                vsr(Colf, Gtc=diag(2)), # diagonal structure at trait level
              rcov = ~ vsr(units, Gtc = unsm(2)), # unstructured at trait level
              data=DT)
#> iteration    LogLik     wall    cpu(sec)   restrained
#>     1      -375.872   21:17:55      1           0
#>     2      -291.932   21:17:56      2           0
#>     3      -258.273   21:17:58      4           0
#>     4      -253.459   21:17:59      5           0
#>     5      -253.291   21:18:0      6           0
#>     6      -253.278   21:18:1      7           0
#>     7      -253.277   21:18:3      9           0
#>     8      -253.277   21:18:4      10           0
summary(mix2)
#> $groups
#>        Yield color
#> u:id     363   363
#> u:Rowf    13    13
#> u:Colf    36    36
#> 
#> $varcomp
#>                          VarComp    VarCompSE     Zratio Constraint
#> u:id.Yield-Yield    7.831409e+02 3.197432e+02  2.4492810   Positive
#> u:id.Yield-color    2.423506e-01 4.134936e-01  0.5861048   Unconstr
#> u:id.color-color    4.921614e-03 9.890881e-04  4.9759104   Positive
#> u:Rowf.Yield-Yield  8.418103e+02 3.939004e+02  2.1371147   Positive
#> u:Rowf.color-color  1.751717e-04 1.243843e-04  1.4083108   Positive
#> u:Colf.Yield-Yield  1.686210e+02 1.232510e+02  1.3681109   Positive
#> u:Colf.color-color  1.844771e-05 9.154208e-05  0.2015216   Positive
#> u:units.Yield-Yield 3.019208e+03 2.837930e+02 10.6387695   Positive
#> u:units.Yield-color 3.649942e-01 1.983709e-01  1.8399580   Unconstr
#> u:units.color-color 2.444059e-03 2.903263e-04  8.4183185   Positive
#> 
#> $betas
#>   Trait      Effect    Estimate   Std.Error  t.value
#> 1 Yield (Intercept) 132.3450577 8.836504001 14.97708
#> 2 color (Intercept)   0.1821964 0.004566864 39.89530
#> 
#> $method
#> [1] "NR"
#> 
#> $logo
#>          logLik      AIC      BIC Method Converge
#> Value -253.2767 510.5535 519.7175     NR     TRUE
#> 
#> attr(,"class")
#> [1] "summary.mmer" "list"        

# }

```
