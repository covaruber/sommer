# **So**lving **M**ixed **M**odel **E**quations in **R** ![Figure: mai.png](figures/mai.png)

Sommer is a structural univariate and multivariate linear mixed-model
package for fitting models with multiple random effects and flexible
covariance structures. Variance parameters are estimated by restricted
maximum likelihood (REML). The package provides two complementary
numerical formulations.

[`mmes`](https://covaruber.github.io/sommer/reference/mmes.md) uses
Henderson's mixed model equations with an Average Information REML
algorithm. The mixed-model coefficient matrix is represented as a sparse
system and factorized with Eigen/Armadillo routines. Sparse LDLT
factorizations are reused throughout the REML calculations. Selected
inverse elements required for likelihood derivatives and prediction
error variances can be obtained with Takahashi sparse-inverse
recursions, so a complete dense coefficient-matrix inverse is not
routinely formed.

[`mmer`](https://covaruber.github.io/sommer/reference/mmer.md) provides
the marginal-covariance MNR formulation. Its REML calculations are
organized around the covariance matrix of the observations and the
corresponding REML projection matrix, with Newton-Raphson/Fisher-scoring
and Average Information updates for variance parameters.

Both formulations support multiple random effects and structured
covariance models. Through
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md), covariance
structures can be composed from an arbitrary number of covariance
factors, including identity, diagonal, unstructured, autoregressive,
moving-average, Toeplitz, factor-analytic, reduced-rank, antedependence,
correlation and spatial structures, as well as user-defined covariance
functions.

The marginal formulation in `mmer` performs its principal REML
calculations in observation space, whereas `mmes` performs them through
the sparse mixed-model coefficient system. Their relative computational
performance therefore depends on the number of observations, the number
of mixed-model coefficients, sparsity, and the covariance structures in
the fitted model.

The numerical algorithms are coded primarily in C++ using Armadillo and
Eigen. Sommer returns REML variance-covariance estimates, BLUEs, BLUPs,
residuals, fitted values, information matrices and, when requested,
prediction error variance information.

## Author

Giovanny Covarrubias-Pazaran

## Functions for genetic analysis

The package provides kernels to estimate additive
([`A.mat`](https://covaruber.github.io/sommer/reference/A.mat.md)),
dominance
([`D.mat`](https://covaruber.github.io/sommer/reference/D.mat.md)),
epistatic
([`E.mat`](https://covaruber.github.io/sommer/reference/E.mat.md)),
single-step
([`H.mat`](https://covaruber.github.io/sommer/reference/H.mat.md))
relationship matrices for diploid and polyploid organisms. It also
provides flexibility to fit other genetic models such as full and half
diallel models and random regression models.

A good converter from letter code to numeric format is implemented in
the function
[`atcg1234`](https://rdrr.io/pkg/enhancer/man/atcg1234.html), which
supports higher ploidy levels than diploid. Additional functions for
genetic analysis have been included such as build a genotypic hybrid
marker matrix
([`build.HMM`](https://rdrr.io/pkg/enhancer/man/build.HMM.html)), plot
of genetic maps
([`map.plot`](https://rdrr.io/pkg/enhancer/man/map.plot.html)), creation
of manhattan plots
([`manhattan`](https://rdrr.io/pkg/enhancer/man/manhattan.html)). If you
need to use pedigree you need to convert your pedigree into a
relationship matrix (use the \`getA\` function from the pedigreemm
package).

## Functions for statistical analysis and S3 methods

The
[`vpredict`](https://covaruber.github.io/sommer/reference/vpredict.md)
function can be used to estimate standard errors for linear combinations
of variance components (e.g. ratios like h2). The
[`r2`](https://covaruber.github.io/sommer/reference/r2.md) function
calculates reliability. S3 methods are available for some parameter
extraction such as:

\+
[`predict.mmes`](https://covaruber.github.io/sommer/reference/predict_mmes.md)

\+
[`fitted.mmes`](https://covaruber.github.io/sommer/reference/fitted_mmes.md)

\+
[`residuals.mmes`](https://covaruber.github.io/sommer/reference/residuals_mmes.md)

\+
[`summary.mmes`](https://covaruber.github.io/sommer/reference/summary_mmes.md)

\+
[`coef.mmes`](https://covaruber.github.io/sommer/reference/coef_mmes.md)

\+
[`anova.mmes`](https://covaruber.github.io/sommer/reference/anova_mmes.md)

\+
[`plot.mmes`](https://covaruber.github.io/sommer/reference/plot_mmes.md)

## Functions for trial analysis

Recently, spatial modeling has been added added to sommer using the
two-dimensional spline
([`spl2Dc`](https://covaruber.github.io/sommer/reference/spl2Dc.md)).

## Keeping sommer updated

The sommer package is updated on CRAN every 4-months due to CRAN
policies but you can find the latest source at
https://github.com/covaruber/sommer. This can be easily installed typing
the following in the R console:

[`library(devtools)`](https://devtools.r-lib.org/)

`install_github("covaruber/sommer")`

This is recommended if you reported a bug, was fixed and was immediately
pushed to GitHub but not in CRAN until the next update.

## Tutorials

**For tutorials** on how to perform different analysis with sommer
please look at the vignettes by typing in the terminal:

[`vignette("sommer.qg")`](https://covaruber.github.io/sommer/articles/sommer.qg.md)

[`vignette("sommer.gxe")`](https://covaruber.github.io/sommer/articles/sommer.gxe.md)

[`vignette("sommer.vs.lme4")`](https://covaruber.github.io/sommer/articles/sommer.vs.lme4.md)

[`vignette("sommer.spatial")`](https://covaruber.github.io/sommer/articles/sommer.spatial.md)

## Getting started

The package has been equiped with several datasets to learn how to use
the sommer package (and almost to learn all sort of quantitative genetic
analysis):

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

## Differences of sommer \>= 4.4.1 with previous versions

Since version 4.4.1, I have unified the use of the two different solving
algorithms into the mmes function by just using the new argument
`henderson` which by default is set to FALSE. Other than that the rest
is the same with the addition that now the identity terms needs to be
encapsulated in the
[`ism`](https://covaruber.github.io/sommer/reference/ism.md) function.
In addition, now the multi-trait models need to be fitted in the long
format. This are few but major changes to the way sommer models are
fitted.

## Differences of sommer \>= 4.1.7 with previous versions

Since version 4.1.7 I have introduced the mmes-based average information
function \`mmec\` which is much faster when dealing with the r \> c
problem (more records than coefficients to estimate). This introduces
its own covariance structure functons such as vsc(), usc(), dsc(),
atc(), csc(). Please give it a try, although is in early phase of
development.

## Differences of sommer \>= 3.7.0 with previous versions

Since version 3.7 I have completly redefined the specification of the
variance-covariance structures to provide more flexibility to the user.
This has particularly helped the residual covariance structures and the
easier combination of custom random effects and overlay models. I think
that although this will bring some uncomfortable situations at the
beggining, in the long term this will help users to fit better models.
In esence, I have abandoned the asreml formulation (not the structures
available) given it's limitations to combine some of the sommer
structures but all covariance structures can now be fitted using the
\`vsm\` functions.

## Differences of sommer \>= 3.0.0 with previous versions

Since version 3.0 I have decided to focus in developing the multivariate
solver and for doing this I have decided to remove the M argument (for
GWAS analysis) from the mmes function and move it to it's own function
GWAS.

Before the mmes solver had implemented the usm(trait), diag(trait),
at(trait) asreml formulation for multivariate models that allow to
specify the structure of the trait in multivariate models. Therefore the
MVM argument was no longer needed. After version 3.7 now the multi-trait
structures can be specified in the `Gt` and `Gtc` arguments of the
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md) function.

The Average Information algorithm had been removed in the past from the
package because of its instability to deal with very complex models
without good initial values. Now after 3.7 I have brought it back after
I noticed that starting with NR the first three iterations gives enough
flexibility to the AI algorithm.

Keep in mind that sommer uses direct inversion (DI) algorithm which can
be very slow for datasets with many observations (big 'n'). The package
is focused in problems of the type p \> n (more random effect(s) levels
than observations) and models with dense covariance structures. For
example, for experiment with dense covariance structures with
low-replication (i.e. 2000 records from 1000 individuals replicated
twice with a covariance structure of 1000x1000) sommer will be faster
than MME-based software. Also for genomic problems with large number of
random effect levels, i.e. 300 individuals (n) with 100,000 genetic
markers (p). On the other hand, for highly replicated trials with small
covariance structures or n \> p (i.e. 2000 records from 200 individuals
replicated 10 times with covariance structure of 200x200) asreml or
other MME-based algorithms will be much faster and I recommend you to
use that software.

## Models Enabled

**General linear mixed model**

Both numerical engines fit the Gaussian linear mixed-model family

\$\$y = X\beta + Zu + e,\$\$

with

\$\$u \sim N(0,G), \qquad e \sim N(0,R),\$\$

and therefore

\$\$V = \mathrm{Var}(y) = ZGZ^\prime + R.\$\$

For several independent random terms this becomes

\$\$V = \sum\_{k=1}^{q} Z_k G_k Z_k^\prime + R.\$\$

Here \\y\\ is the response vector, \\X\\ and \\Z\\ are design matrices,
\\\beta\\ contains fixed effects, \\u\\ contains random effects, and
\\e\\ contains residual effects. The same representation covers
univariate and multivariate models after the corresponding responses,
design matrices and covariance structures are assembled.

Covariance models specified through
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md) use a
product-level variance scale and an arbitrary number of
covariance-shaping factors:

\$\$\Sigma = \sigma^2(K_1 \otimes K_2 \otimes \cdots \otimes K_m).\$\$

For a random effect with known relationship covariance \\A\\,

\$\$G = \sigma^2(K_1 \otimes K_2 \otimes \cdots \otimes K_m) \otimes
A.\$\$

The CovarianceFactor interface separates evaluation of each \\K_j\\, its
derivatives, parameter reporting and trust limits from the REML solver.
New covariance structures can therefore be added without changing the
core `mmes` optimization algorithm.

**Marginal MNR formulation in mmer**

The MNR engine used by
[`mmer`](https://covaruber.github.io/sommer/reference/mmer.md) organizes
the REML calculations around the marginal covariance \\V\\. Define the
REML projection matrix

\$\$P = V^{-1} - V^{-1}X(X^\prime V^{-1}X)^{-1}X^\prime V^{-1}.\$\$

Apart from constants independent of the variance parameters, the
Gaussian REML log-likelihood is

\$\$\ell_R(\theta) = -\frac{1}{2}\\\log\|V\|+\log\|X^\prime
V^{-1}X\|+y^\prime P y\\.\$\$

For variance parameter \\\theta_i\\, let

\$\$V_i = \frac{\partial V}{\partial\theta_i}.\$\$

For covariance models linear in the fitted variance parameters, the REML
score has the standard form

\$\$s_i = \frac{1}{2}\\y^\prime P V_i P y-\mathrm{tr}(P V_i)\\.\$\$

The MNR C++ implementation forms the projected covariance derivatives
\\PV_i\\. Its Average Information matrix is

\$\$AI\_{ij} = \frac{1}{2}y^\prime P V_i P V_j P y,\$\$

while its Newton-Raphson/Fisher-scoring path uses trace products of the
projected covariance derivatives. The resulting information system
determines the variance-parameter search direction, and the
implementation can blend the information calculation with an EM-type
diagonal stabilization during early iterations.

This formulation is referred to as *marginal* because its principal
covariance calculations are carried out in observation space through
\\V\\ and \\P\\. This describes the numerical formulation more
accurately than naming it after a particular matrix operation.

**Sparse Henderson Average Information formulation in mmes**

The [`mmes`](https://covaruber.github.io/sommer/reference/mmes.md)
engine uses Henderson's mixed model equations. With \\W=\[X\\Z\]\\, the
coefficient system is

\$\$C = \left\[ \begin{array}{cc} X^\prime R^{-1}X & X^\prime R^{-1}Z \\
Z^\prime R^{-1}X & Z^\prime R^{-1}Z+G^{-1} \end{array} \right\],\$\$

and the right-hand side is

\$\$r = \left\[ \begin{array}{c} X^\prime R^{-1}y \\ Z^\prime R^{-1}y
\end{array} \right\].\$\$

The BLUE and BLUP solutions satisfy

\$\$Cb=r, \qquad
b=(\widehat{\beta}^\prime,\widehat{u}^\prime)^\prime.\$\$

The current C++ implementation stores \\C\\ sparsely and uses an Eigen
`SimplicialLDLT` factorization. Symbolic factorization information is
reused when the sparsity pattern remains unchanged. Random-effect
precision contributions are inserted blockwise. Structured residual
covariance models use optimized diagonal or repeated-block/Kronecker
paths when applicable and a generic sparse residual-factorization path
otherwise.

The REML criterion is evaluated from determinant contributions for the
random covariance structures, residual covariance and mixed-model
coefficient matrix, together with the quadratic term \\y^\prime P y\\.
This is a mixed-model-equation representation of the same REML objective
optimized by the marginal formulation.

For covariance parameter \\\eta_i\\, the generic covariance engine
supplies

\$\$\Sigma_i = \frac{\partial\Sigma}{\partial\eta_i}.\$\$

For

\$\$\Sigma=\sigma^2(K_1\otimes\cdots\otimes K_m),\$\$

a parameter belonging to factor \\K_j\\ has derivative

\$\$\frac{\partial\Sigma}{\partial\eta\_{jk}} = \sigma^2
K_1\otimes\cdots\otimes \frac{\partial K_j}{\partial\eta\_{jk}}
\otimes\cdots\otimes K_m.\$\$

For the working scale \\\tau=\log(\sigma^2)\\,

\$\$\frac{\partial\Sigma}{\partial\tau}=\Sigma.\$\$

The score trace terms use selected elements associated with the sparse
coefficient system. These are obtained from the LDLT factors with
Takahashi sparse-inverse recursions; a complete \\C^{-1}\\ is not
required during REML iterations.

Average Information is calculated by analytically differentiating the
factorized mixed-model equations. Differentiating

\$\$Cb=r\$\$

with respect to \\\eta_i\\ gives

\$\$C b_i = r_i-C_i b,\$\$

where

\$\$C_i=\frac{\partial C}{\partial\eta_i}, \qquad r_i=\frac{\partial
r}{\partial\eta_i}.\$\$

The sensitivity systems therefore reuse the same sparse factorization of
\\C\\. For covariance \\\Sigma\\ with precision \\\Lambda=\Sigma^{-1}\\
and derivative \\B_i=\partial\Sigma/\partial\eta_i\\,

\$\$\Lambda_i=-\Lambda B_i\Lambda.\$\$

This factor-differentiation formulation avoids constructing an
observation-space working variate for every random covariance parameter.
Residual/residual and random/residual Average Information terms use the
corresponding mixed-model/Schur-complement identities.

Covariance constructors may provide analytic derivatives in C++ or R. If
a constructor uses numerical differentiation, only its small covariance
factor \\K_j(\eta)\\ is differentiated by central finite differences;
the REML likelihood and downstream score and Average Information
calculations are not numerically differentiated.

**Sparse inverse subsets and prediction error variances**

During optimization, `mmes` does not routinely materialize the complete
\\C^{-1}\\. Takahashi recursions provide the selected inverse subset
required for REML trace calculations. After convergence, `computeCi=1`
obtains the diagonal information needed for prediction error variances
through the sparse inverse subset, whereas `computeCi=2` requests the
complete coefficient-matrix inverse.

**Parameter updates and numerical safeguards**

Covariance parameters are optimized in working coordinates chosen to
respect their parameter spaces, such as logarithms for positive scales
and hyperbolic-tangent or bounded-logit mappings for correlations. Each
CovarianceFactor descriptor supplies its natural-scale reporting
transformation and parameter-specific trust limit.

The information system determines a proposed REML update. In the current
`mmes` implementation, the proposal is limited in working-coordinate
space and then checked against the exact REML likelihood. If the
proposed global step decreases the accepted likelihood beyond numerical
tolerance, the step is geometrically reduced and reevaluated. A
likelihood decrease is therefore not interpreted as convergence.
Boundary handling and robust information solves are retained for
difficult or weakly identified variance components.

**Choosing between the two formulations**

Both engines estimate the same class of linear mixed models but organize
the numerical work differently. `mmer` uses the marginal covariance
representation and is attractive when observation-space calculations are
moderate in size. `mmes` uses the sparse Henderson coefficient system
and is particularly attractive when the mixed-model equations are sparse
and their factorization can be reused efficiently. The best choice
depends on the dimensions, sparsity and covariance structures of the
fitted model rather than on a universal rule based only on the number of
records or coefficients.

See [`mmes`](https://covaruber.github.io/sommer/reference/mmes.md),
[`mmer`](https://covaruber.github.io/sommer/reference/mmer.md), and
[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md) for the
interfaces and available covariance structures.

## Bug report and contact

If you have any questions or suggestions please post it in
https://stackoverflow.com or https://stats.stackexchange.com

I'll be glad to help or answer any question. I have spent a valuable
amount of time developing this package. Please cite this package in your
publication. Type 'citation("sommer")' to know how to cite it.

## References

Covarrubias-Pazaran G. 2016. Genome assisted prediction of quantitative
traits using the R package sommer. PLoS ONE 11(6):
doi:10.1371/journal.pone.0156744

Covarrubias-Pazaran G. 2018. Software update: Moving the R package
sommer to multivariate mixed models for genome-assisted prediction. doi:
https://doi.org/10.1101/354639

Sanderson, C., & Curtin, R. (2025). Armadillo: An Efficient Framework
for Numerical Linear Algebra. arXiv preprint arXiv:2502.03000.

Bernardo Rex. 2010. Breeding for quantitative traits in plants. Second
edition. Stemma Press. 390 pp.

Gilmour et al. 1995. Average Information REML: An efficient algorithm
for variance parameter estimation in linear mixed models. Biometrics
51(4):1440-1450.

Henderson C.R. 1975. Best Linear Unbiased Estimation and Prediction
under a Selection Model. Biometrics vol. 31(2):423-447.

Kang et al. 2008. Efficient control of population structure in model
organism association mapping. Genetics 178:1709-1723.

Lee et al. 2015. MTG2: An efficient algorithm for multivariate linear
mixed model analysis based on genomic information. Cold Spring Harbor.
doi: http://dx.doi.org/10.1101/027201.

Maier et al. 2015. Joint analysis of psychiatric disorders increases
accuracy of risk prediction for schizophrenia, bipolar disorder, and
major depressive disorder. Am J Hum Genet; 96(2):283-294.

Searle. 1993. Applying the EM algorithm to calculating ML and REML
estimates of variance components. Paper invited for the 1993 American
Statistical Association Meeting, San Francisco.

Yu et al. 2006. A unified mixed-model method for association mapping
that accounts for multiple levels of relatedness. Genetics 38:203-208.

Tunnicliffe W. 1989. On the use of marginal likelihood in time series
model estimation. JRSS 51(1):15-27.

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

# \donttest{



DT <- DT_example
DT=DT[with(DT, order(Env)), ]
head(DT)
#>          Name     Env Loc Year     Block Yield      Weight
#> 67   MSL007-B CA.2011  CA 2011 CA.2011.2     5 -1.43551031
#> 105  MSL007-B CA.2011  CA 2011 CA.2011.1     6 -1.43896249
#> 308  MSK061-4 CA.2011  CA 2011 CA.2011.2     9 -0.79839374
#> 393  MSK061-4 CA.2011  CA 2011 CA.2011.1    10 -0.19022655
#> 469 MSR169-8Y CA.2011  CA 2011 CA.2011.1    11  0.02682037
#> 471     NY148 CA.2011  CA 2011 CA.2011.1    11 -0.03278206

####=========================================####
#### Univariate homogeneous variance models  ####
####=========================================####

## Compound simmetry (CS) model
ans1 <- mmes(Yield~Env,
             random= ~ Name + Env:Name,
             rcov= ~ units,
             data=DT)
#> Solver selected: ldlt
#> OpenMP available: up to 8 threads.
#> OpenMP active: parallel SLQ log-determinant probes (8 probes).
#> iteration    LogLik     wall    cpu(sec)   restrained   EM weight      pivot
#>     1      -31.6987   21:18:7      0           0      1      3.08144
#>     2      -28.0631   21:18:7      0           0      0.813615      5.42039
#>     3      -26.3167   21:18:7      0           0      0.661969      3.92118
#>     4      -26.2862   21:18:7      0           0      0.538588      4.00163
#>     5      -26.2851   21:18:7      0           0      0.438203      3.9582
#>     6      -26.2851   21:18:7      0           0      0.356529      3.95514
summary(ans1)
#> ============================================================
#>          Multivariate Linear Mixed Model fit by  REML         
#> **********************  sommer 4.4  ********************** 
#> ============================================================
#>          logLik      AIC      BIC Method Converge
#> Value -26.28509 58.57018 68.23125     AI     TRUE
#> ============================================================
#> Variance-Covariance components:
#>                 term factor parameter estimate StdError Zratio
#> 1     vsm(ism(Name)) sigma2    sigma2    3.688   1.1080  3.329
#> 2 vsm(ism(Env:Name)) sigma2    sigma2    5.167   0.9878  5.230
#> 3    vsm(ism(units)) sigma2    sigma2    4.369   0.4508  9.692
#> ============================================================
#> Fixed effects:
#>            Estimate Std.Error t.value
#> Intercept    16.496        NA      NA
#> EnvCA.2012   -5.776        NA      NA
#> EnvCA.2013   -6.380        NA      NA
#> ============================================================
#> Use the '$' sign to access results and parameters

####===========================================####
#### Univariate heterogeneous variance models  ####
####===========================================####
## Compound simmetry (CS) + Diagonal (DIAG) model
ans3 <- mmes(Yield~Env,
             random= ~Name + vsm(dsm(Env),ism(Name)),
             rcov= ~ vsm(dsm(Env),ism(units)),
             data=DT)
#> Solver selected: ldlt
#> OpenMP available: up to 8 threads.
#> OpenMP active: parallel SLQ log-determinant probes (8 probes).
#>   Working-coordinate trust scaling applied: alpha=0.650987
#> iteration    LogLik     wall    cpu(sec)   restrained   EM weight      pivot
#>     1      -31.6987   21:18:7      0           0      1      3.08144
#>     2      -29.0398   21:18:7      0           0      0.813615      4.45068
#>     3      -24.0636   21:18:7      0           0      0.661969      2.4284
#>     4      -22.7208   21:18:7      0           0      0.538588      3.22302
#>     5      -21.8476   21:18:7      0           0      0.438203      2.6077
#>     6      -21.7347   21:18:7      0           0      0.356529      2.54008
#>     7      -21.667   21:18:7      0           0      0.290077      2.4098
#>     8      -21.616   21:18:7      0           0      0.236011      2.2611
#>     9      -21.5897   21:18:7      0           0      0.192022      2.15979
#>     10      -21.5774   21:18:7      0           0      0.156232      2.09536
#>     11      -21.5723   21:18:7      0           0      0.127113      2.05671
#>     12      -21.5703   21:18:7      0           0      0.103421      2.03495
#>     13      -21.5697   21:18:7      0           0      0.0841447      2.02349
#>     14      -21.5696   21:18:7      0           0      0.0684614      2.01786
#>     15      -21.5695   21:18:7      0           0      0.0557012      2.01529
summary(ans3)
#> ============================================================
#>          Multivariate Linear Mixed Model fit by  REML         
#> **********************  sommer 4.4  ********************** 
#> ============================================================
#>          logLik      AIC      BIC Method Converge
#> Value -21.56953 49.13905 58.80012     AI     TRUE
#> ============================================================
#> Variance-Covariance components:
#>                        term factor         parameter estimate StdError Zratio
#> 1            vsm(ism(Name)) sigma2            sigma2    2.963   1.0153  2.918
#> 2  vsm(dsm(Env), ism(Name))   diag variance[CA.2011]   10.140   2.9349  3.455
#> 3  vsm(dsm(Env), ism(Name))   diag variance[CA.2012]    1.879   1.3123  1.432
#> 4  vsm(dsm(Env), ism(Name))   diag variance[CA.2013]    6.630   1.7626  3.762
#> 5 vsm(dsm(Env), ism(units))   diag variance[CA.2011]    4.948   0.9092  5.442
#> 6 vsm(dsm(Env), ism(units))   diag variance[CA.2012]    5.724   0.9345  6.125
#> 7 vsm(dsm(Env), ism(units))   diag variance[CA.2013]    2.559   0.4525  5.656
#> ============================================================
#> Fixed effects:
#>            Estimate Std.Error t.value
#> Intercept    16.508        NA      NA
#> EnvCA.2012   -5.817        NA      NA
#> EnvCA.2013   -6.412        NA      NA
#> ============================================================
#> Use the '$' sign to access results and parameters



# }
```
