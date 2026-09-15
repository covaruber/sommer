\name{sommer-package}
\alias{sommer}
\docType{package}
\author{
Giovanny Covarrubias-Pazaran
}
\title{

\strong{So}lving \strong{M}ixed \strong{M}odel \strong{E}quations in \strong{R}
\cr

\if{html}{\figure{mai.png}{options: style="height: 10\%;" alt="Figure: mai.png"}}

}
\description{

Sommer is a structural univariate and multivariate linear mixed-model package
for fitting models with multiple random effects and flexible covariance
structures. Variance parameters are estimated by restricted maximum likelihood
(REML). The package provides two complementary numerical formulations.

\code{\link{mmes}} uses Henderson's mixed model equations with an Average
Information REML algorithm. The mixed-model coefficient matrix is represented
as a sparse system and factorized with Eigen/Armadillo routines. Sparse LDLT
factorizations are reused throughout the REML calculations. Selected inverse
elements required for likelihood derivatives and prediction error variances can
be obtained with Takahashi sparse-inverse recursions, so a complete dense
coefficient-matrix inverse is not routinely formed.

\code{\link{mmer}} provides the marginal-covariance MNR formulation. Its REML
calculations are organized around the covariance matrix of the observations and
the corresponding REML projection matrix, with Newton-Raphson/Fisher-scoring
and Average Information updates for variance parameters.

Both formulations support multiple random effects and structured covariance
models. Through \code{\link{vsm}}, covariance structures can be composed from
an arbitrary number of covariance factors, including identity, diagonal,
unstructured, autoregressive, moving-average, Toeplitz, factor-analytic,
reduced-rank, antedependence, correlation and spatial structures, as well as
user-defined covariance functions.

The marginal formulation in \code{mmer} performs its principal REML
calculations in observation space, whereas \code{mmes} performs them through
the sparse mixed-model coefficient system. Their relative computational
performance therefore depends on the number of observations, the number of
mixed-model coefficients, sparsity, and the covariance structures in the
fitted model.

The numerical algorithms are coded primarily in C++ using Armadillo and Eigen.
Sommer returns REML variance-covariance estimates, BLUEs, BLUPs, residuals,
fitted values, information matrices and, when requested, prediction error
variance information.

}

\section{Functions for genetic analysis}{

The package provides kernels to estimate additive (\code{\link{A.mat}}), dominance
(\code{\link{D.mat}}), epistatic (\code{\link{E.mat}}), single-step (\code{\link{H.mat}}) 
relationship matrices for diploid and polyploid organisms. It also provides flexibility 
to fit other genetic models such as full and half diallel models and random regression models.

A good converter from letter code to numeric format is implemented in the function 
\code{\link[enhancer]{atcg1234}}, which supports higher ploidy levels than diploid. 
Additional functions for genetic analysis have been included such as build a 
genotypic hybrid marker matrix (\code{\link[enhancer]{build.HMM}}), plot of genetic maps
(\code{\link[enhancer]{map.plot}}), creation of manhattan plots (\code{\link[enhancer]{manhattan}}).
If you need to use pedigree you need to convert your pedigree into a relationship 
matrix (use the `getA` function from the pedigreemm package).

}

\section{Functions for statistical analysis and S3 methods}{

The \code{\link{vpredict}} function can be used to estimate standard errors for linear combinations of variance components (e.g. ratios like h2). The \code{\link{r2}} function calculates reliability. S3 methods are available for some parameter extraction such as:

+ \code{\link{predict.mmes}}

+ \code{\link{fitted.mmes}}

+ \code{\link{residuals.mmes}}

+ \code{\link{summary.mmes}}

+ \code{\link{coef.mmes}}

+ \code{\link{anova.mmes}}

+ \code{\link{plot.mmes}}

}

\section{Functions for trial analysis}{

Recently, spatial modeling has been added added to sommer using the two-dimensional 
spline (\code{\link{spl2Dc}}).

}

\section{Keeping sommer updated}{
The sommer package is updated on CRAN every 4-months due to CRAN policies but you can find the latest source at https://github.com/covaruber/sommer. This can be easily installed typing the following in the R console:

\code{library(devtools)}

\code{install_github("covaruber/sommer")}

This is recommended if you reported a bug, was fixed and was immediately pushed to GitHub but not in CRAN until the next update.
}

\section{Tutorials}{

 \strong{For tutorials} on how to perform different analysis with sommer please look at the vignettes by typing in the terminal:

\code{vignette("sommer.qg")}

\code{vignette("sommer.gxe")}

\code{vignette("sommer.vs.lme4")} 

\code{vignette("sommer.spatial")}

}

\section{Getting started}{
The package has been equiped with several datasets to learn how to use the sommer package (and almost to learn all sort of quantitative genetic analysis): 

* \code{\link[enhancer]{DT_halfdiallel}}, \code{\link[enhancer]{DT_fulldiallel}} and \code{\link[enhancer]{DT_mohring}} datasets have examples to fit half and full diallel designs. 

* \code{\link[enhancer]{DT_h2}} to calculate heritability

* \code{\link[enhancer]{DT_cornhybrids}} and \code{\link[enhancer]{DT_technow}} datasets to perform genomic prediction in hybrid single crosses 

* \code{\link[enhancer]{DT_wheat}} dataset to do genomic prediction in single crosses in species displaying only additive effects.

* \code{\link[enhancer]{DT_cpdata}} dataset to fit genomic prediction models within a biparental population coming from 2 highly heterozygous parents including additive, dominance and epistatic effects. 

* \code{\link[enhancer]{DT_polyploid}} to fit genomic prediction and GWAS analysis in polyploids. 

* \code{\link[enhancer]{DT_gryphon}} data contains an example of an animal model including pedigree information.

* \code{\link[enhancer]{DT_btdata}} dataset contains an animal (birds) model.

* \code{\link[enhancer]{DT_legendre}} simulated dataset for random regression model.

* \code{\link[enhancer]{DT_sleepstudy}} dataset to know how to translate lme4 models to sommer models.

}

\section{Differences of sommer >= 4.4.1 with previous versions}{

Since version 4.4.1, I have unified the use of the two different solving algorithms
into the mmes function by just using the new argument \code{henderson} which by default 
is set to FALSE. Other than that the rest is the same with the addition that now
the identity terms needs to be encapsulated in the \code{\link{ism}} function. In 
addition, now the multi-trait models need to be fitted in the long format. This are 
few but major changes to the way sommer models are fitted.

}

\section{Differences of sommer >= 4.1.7 with previous versions}{

Since version 4.1.7 I have introduced the mmes-based average information function `mmec` which is much faster when dealing with the r > c problem (more records than coefficients to estimate). This introduces its own covariance structure functons such as vsc(), usc(), dsc(), atc(), csc(). Please give it a try, although is in early phase of development.

}

\section{Differences of sommer >= 3.7.0 with previous versions}{

Since version 3.7 I have completly redefined the specification of the variance-covariance structures to provide more flexibility to the user. This has particularly helped the residual covariance structures and the easier combination of custom random effects and overlay models. I think that although this will bring some uncomfortable situations at the beggining, in the long term this will help users to fit better models. In esence, I have abandoned the asreml formulation (not the structures available) given it's limitations to combine some of the sommer structures but all covariance structures can now be fitted using the `vsm` functions. 

}
\section{Differences of sommer >= 3.0.0 with previous versions}{

Since version 3.0 I have decided to focus in developing the multivariate solver and for doing this I have decided to remove the M argument (for GWAS analysis) from the mmes function and move it to it's own function GWAS.

Before the mmes solver had implemented the usm(trait), diag(trait), at(trait) asreml formulation for multivariate models that allow to specify the structure of the trait in multivariate models. Therefore the MVM argument was no longer needed. After version 3.7 now the multi-trait structures can be specified in the \code{Gt} and \code{Gtc} arguments of the \code{\link{vsm}} function.

The Average Information algorithm had been removed in the past from the package because of its instability to deal with very complex models without good initial values. Now after 3.7 I have brought it back after I noticed that starting with NR the first three iterations gives enough flexibility to the AI algorithm.

Keep in mind that sommer uses direct inversion (DI) algorithm which can be very slow for datasets with many observations (big 'n'). The package is focused in problems of the type p > n (more random effect(s) levels than observations) and models with dense covariance structures. For example, for experiment with dense covariance structures with low-replication (i.e. 2000 records from 1000 individuals replicated  twice with a covariance structure of 1000x1000) sommer will be faster than MME-based software. Also for genomic problems with large number of random effect levels, i.e. 300 individuals (n) with 100,000 genetic markers (p). On the other hand, for highly replicated trials with small covariance structures or n > p (i.e. 2000 records from 200 individuals replicated 10 times with covariance structure of 200x200) asreml or other MME-based algorithms will be much faster and I recommend you to use that software.

}

\section{Models Enabled}{

\strong{General linear mixed model}

Both numerical engines fit the Gaussian linear mixed-model family

\deqn{y = X\beta + Zu + e,}

with

\deqn{u \sim N(0,G), \qquad e \sim N(0,R),}

and therefore

\deqn{V = \mathrm{Var}(y) = ZGZ^\prime + R.}

For several independent random terms this becomes

\deqn{V = \sum_{k=1}^{q} Z_k G_k Z_k^\prime + R.}

Here \eqn{y} is the response vector, \eqn{X} and \eqn{Z} are design matrices,
\eqn{\beta} contains fixed effects, \eqn{u} contains random effects, and
\eqn{e} contains residual effects. The same representation covers univariate
and multivariate models after the corresponding responses, design matrices and
covariance structures are assembled.

Covariance models specified through \code{\link{vsm}} use a product-level
variance scale and an arbitrary number of covariance-shaping factors:

\deqn{\Sigma = \sigma^2(K_1 \otimes K_2 \otimes \cdots \otimes K_m).}

For a random effect with known relationship covariance \eqn{A},

\deqn{G = \sigma^2(K_1 \otimes K_2 \otimes \cdots \otimes K_m) \otimes A.}

The CovarianceFactor interface separates evaluation of each \eqn{K_j}, its
derivatives, parameter reporting and trust limits from the REML solver. New
covariance structures can therefore be added without changing the core
\code{mmes} optimization algorithm.

\strong{Marginal MNR formulation in mmer}

The MNR engine used by \code{\link{mmer}} organizes the REML calculations
around the marginal covariance \eqn{V}. Define the REML projection matrix

\deqn{P = V^{-1} -
V^{-1}X(X^\prime V^{-1}X)^{-1}X^\prime V^{-1}.}

Apart from constants independent of the variance parameters, the Gaussian REML
log-likelihood is

\deqn{\ell_R(\theta) =
-\frac{1}{2}\{\log|V|+\log|X^\prime V^{-1}X|+y^\prime P y\}.}

For variance parameter \eqn{\theta_i}, let

\deqn{V_i = \frac{\partial V}{\partial\theta_i}.}

For covariance models linear in the fitted variance parameters, the REML score
has the standard form

\deqn{s_i =
\frac{1}{2}\{y^\prime P V_i P y-\mathrm{tr}(P V_i)\}.}

The MNR C++ implementation forms the projected covariance derivatives
\eqn{PV_i}. Its Average Information matrix is

\deqn{AI_{ij} =
\frac{1}{2}y^\prime P V_i P V_j P y,}

while its Newton-Raphson/Fisher-scoring path uses trace products of the
projected covariance derivatives. The resulting information system determines
the variance-parameter search direction, and the implementation can blend the
information calculation with an EM-type diagonal stabilization during early
iterations.

This formulation is referred to as \emph{marginal} because its principal
covariance calculations are carried out in observation space through \eqn{V}
and \eqn{P}. This describes the numerical formulation more accurately than
naming it after a particular matrix operation.

\strong{Sparse Henderson Average Information formulation in mmes}

The \code{\link{mmes}} engine uses Henderson's mixed model equations. With
\eqn{W=[X\;Z]}, the coefficient system is

\deqn{C =
\left[
\begin{array}{cc}
X^\prime R^{-1}X & X^\prime R^{-1}Z \\
Z^\prime R^{-1}X & Z^\prime R^{-1}Z+G^{-1}
\end{array}
\right],}

and the right-hand side is

\deqn{r =
\left[
\begin{array}{c}
X^\prime R^{-1}y \\
Z^\prime R^{-1}y
\end{array}
\right].}

The BLUE and BLUP solutions satisfy

\deqn{Cb=r, \qquad
b=(\widehat{\beta}^\prime,\widehat{u}^\prime)^\prime.}

The current C++ implementation stores \eqn{C} sparsely and uses an Eigen
\code{SimplicialLDLT} factorization. Symbolic factorization information is
reused when the sparsity pattern remains unchanged. Random-effect precision
contributions are inserted blockwise. Structured residual covariance models
use optimized diagonal or repeated-block/Kronecker paths when applicable and a
generic sparse residual-factorization path otherwise.

The REML criterion is evaluated from determinant contributions for the random
covariance structures, residual covariance and mixed-model coefficient matrix,
together with the quadratic term \eqn{y^\prime P y}. This is a
mixed-model-equation representation of the same REML objective optimized by
the marginal formulation.

For covariance parameter \eqn{\eta_i}, the generic covariance engine supplies

\deqn{\Sigma_i =
\frac{\partial\Sigma}{\partial\eta_i}.}

For

\deqn{\Sigma=\sigma^2(K_1\otimes\cdots\otimes K_m),}

a parameter belonging to factor \eqn{K_j} has derivative

\deqn{\frac{\partial\Sigma}{\partial\eta_{jk}} =
\sigma^2 K_1\otimes\cdots\otimes
\frac{\partial K_j}{\partial\eta_{jk}}
\otimes\cdots\otimes K_m.}

For the working scale \eqn{\tau=\log(\sigma^2)},

\deqn{\frac{\partial\Sigma}{\partial\tau}=\Sigma.}

The score trace terms use selected elements associated with the sparse
coefficient system. These are obtained from the LDLT factors with Takahashi
sparse-inverse recursions; a complete \eqn{C^{-1}} is not required during REML
iterations.

Average Information is calculated by analytically differentiating the
factorized mixed-model equations. Differentiating

\deqn{Cb=r}

with respect to \eqn{\eta_i} gives

\deqn{C b_i = r_i-C_i b,}

where

\deqn{C_i=\frac{\partial C}{\partial\eta_i}, \qquad
r_i=\frac{\partial r}{\partial\eta_i}.}

The sensitivity systems therefore reuse the same sparse factorization of
\eqn{C}. For covariance \eqn{\Sigma} with precision
\eqn{\Lambda=\Sigma^{-1}} and derivative
\eqn{B_i=\partial\Sigma/\partial\eta_i},

\deqn{\Lambda_i=-\Lambda B_i\Lambda.}

This factor-differentiation formulation avoids constructing an
observation-space working variate for every random covariance parameter.
Residual/residual and random/residual Average Information terms use the
corresponding mixed-model/Schur-complement identities.

Covariance constructors may provide analytic derivatives in C++ or R. If a
constructor uses numerical differentiation, only its small covariance factor
\eqn{K_j(\eta)} is differentiated by central finite differences; the REML
likelihood and downstream score and Average Information calculations are not
numerically differentiated.

\strong{Sparse inverse subsets and prediction error variances}

During optimization, \code{mmes} does not routinely materialize the complete
\eqn{C^{-1}}. Takahashi recursions provide the selected inverse subset required
for REML trace calculations. After convergence, \code{computeCi=1} obtains the
diagonal information needed for prediction error variances through the sparse
inverse subset, whereas \code{computeCi=2} requests the complete
coefficient-matrix inverse.

\strong{Parameter updates and numerical safeguards}

Covariance parameters are optimized in working coordinates chosen to respect
their parameter spaces, such as logarithms for positive scales and
hyperbolic-tangent or bounded-logit mappings for correlations. Each
CovarianceFactor descriptor supplies its natural-scale reporting
transformation and parameter-specific trust limit.

The information system determines a proposed REML update. In the current
\code{mmes} implementation, the proposal is limited in working-coordinate
space and then checked against the exact REML likelihood. If the proposed
global step decreases the accepted likelihood beyond numerical tolerance, the
step is geometrically reduced and reevaluated. A likelihood decrease is
therefore not interpreted as convergence. Boundary handling and robust
information solves are retained for difficult or weakly identified variance
components.

\strong{Choosing between the two formulations}

Both engines estimate the same class of linear mixed models but organize the
numerical work differently. \code{mmer} uses the marginal covariance
representation and is attractive when observation-space calculations are
moderate in size. \code{mmes} uses the sparse Henderson coefficient system and
is particularly attractive when the mixed-model equations are sparse and
their factorization can be reused efficiently. The best choice depends on the
dimensions, sparsity and covariance structures of the fitted model rather than
on a universal rule based only on the number of records or coefficients.

See \code{\link{mmes}}, \code{\link{mmer}}, and \code{\link{vsm}} for the
interfaces and available covariance structures.

}

\section{Bug report and contact}{

If you have any questions or suggestions please post it in https://stackoverflow.com or https://stats.stackexchange.com

I'll be glad to help or answer any question. I have spent a valuable amount of time developing this package. Please cite this package in your publication. Type 'citation("sommer")' to know how to cite it.


}
\references{

Covarrubias-Pazaran G. 2016. Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

Covarrubias-Pazaran G. 2018. Software update: Moving the R package sommer to multivariate mixed models for genome-assisted prediction. doi: https://doi.org/10.1101/354639

Sanderson, C., & Curtin, R. (2025). Armadillo: An Efficient Framework for Numerical Linear Algebra. arXiv preprint arXiv:2502.03000.

Bernardo Rex. 2010. Breeding for quantitative traits in plants. Second edition. Stemma Press. 390 pp.

Gilmour et al. 1995. Average Information REML: An efficient algorithm for variance parameter estimation in linear mixed models. Biometrics 51(4):1440-1450.

Henderson C.R. 1975. Best Linear Unbiased Estimation and Prediction under a Selection Model. Biometrics vol. 31(2):423-447.

Kang et al. 2008. Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723.

Lee et al. 2015. MTG2: An efficient algorithm for multivariate linear mixed model analysis based on genomic information. Cold Spring Harbor. doi: http://dx.doi.org/10.1101/027201.

Maier et al. 2015. Joint analysis of psychiatric disorders increases accuracy of risk prediction for schizophrenia, bipolar disorder, and major depressive disorder. Am J Hum Genet; 96(2):283-294.

Searle. 1993. Applying the EM algorithm to calculating ML and REML estimates of variance components. Paper invited for the 1993 American Statistical Association Meeting, San Francisco.

Yu et al. 2006. A unified mixed-model method for association mapping that accounts for multiple levels of relatedness. Genetics 38:203-208.

Tunnicliffe W. 1989. On the use of marginal likelihood in time series model estimation. JRSS 51(1):15-27.

}
\keyword{ R package }
\examples{

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

\donttest{



DT <- DT_example
DT=DT[with(DT, order(Env)), ]
head(DT)

####=========================================####
#### Univariate homogeneous variance models  ####
####=========================================####

## Compound simmetry (CS) model
ans1 <- mmes(Yield~Env,
             random= ~ Name + Env:Name,
             rcov= ~ units,
             data=DT)
summary(ans1)

####===========================================####
#### Univariate heterogeneous variance models  ####
####===========================================####
## Compound simmetry (CS) + Diagonal (DIAG) model
ans3 <- mmes(Yield~Env,
             random= ~Name + vsm(dsm(Env),ism(Name)),
             rcov= ~ vsm(dsm(Env),ism(units)),
             data=DT)
summary(ans3)



}

}