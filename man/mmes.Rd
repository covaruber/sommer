\name{mmes}
\alias{mmes}
\title{

\strong{m}ixed \strong{m}odel \strong{e}quations \strong{s}olver

}
\description{

Fits a generalized linear mixed model
coded in C++ using the Armadillo and Eigen libraries to optimize matrix 
operations. Special covariance structures can be specified with the \code{\link{vsm}()}
function. Two engines are available; henderson formulation and direct reml (henderson argument) and 
three solvers in each; ldlt, cholmod and pcg (solver argument) to adjust to the
modeling situation.


}

\usage{

mmes(fixed, random, rcov, data, W,
     nIters=30, tolParConvLL=1e-04,
     tolParConvNorm=1e-04, tolParInv=1e-06,
     naMethodX="exclude", naMethodY="exclude",
     naMethodRandom="exclude", naMethodR="exclude",
     returnParam=FALSE, dateWarning=TRUE,
     verbose=TRUE, stepWeight=NULL, emWeight=NULL,
     contrasts=NULL, getPEV=TRUE, henderson=TRUE,
    computeCi=0, solver="auto", pcgTol=1.0e-8,
    pcgMaxIters=0, pcgTraceProbes=8,
    pcgLanczosSteps=20, REML=TRUE,
    family=stats::gaussian(), pqlControl=list(),
    .pqlInner=FALSE, .pqlFixedDispersion=FALSE,
    .pqlWorkingPrecision=NULL, .pqlBaseW=NULL,
    .pqlBaseFactor=NULL)
}
\arguments{
  \item{fixed}{A formula specifying the \strong{response variable(s)} \strong{and fixed effects}, i.e:

  \emph{response ~ covariate}

  }

  \item{random}{A formula specifying the \strong{random effects}, e.g.,
  \emph{random = ~ genotype + year}. A simple random term that is not wrapped
  in \code{\link{vsm}} is treated as an identity-covariance random effect.

  The \code{\link{vsm}} function is the main interface for specifying structured
  covariance models for random effects. Its current covariance model is

  \deqn{G = \sigma^2 (K_1 \otimes K_2 \otimes \cdots \otimes K_m) \otimes A,}

  where \eqn{\sigma^2} is the single variance scale owned by \code{vsm()},
  \eqn{K_1,\ldots,K_m} are covariance-shaping factors supplied by covariance
  constructors, and \eqn{A} is the known covariance relationship for the final
  random effect. In the Henderson implementation \code{Gu} supplies the
  \strong{inverse} of \eqn{A}; it must have row and column names matching the
  levels of the final random effect and \code{attr(Gu,"inverse")=TRUE}.

  Within \code{vsm(...)}, the \strong{last} constructor supplies the incidence
  matrix for the main random effect, while all preceding constructors define
  covariance-shaping factors. The number of covariance factors is not limited.
  For example,

  \code{random = ~ vsm(dsm(Location), ism(Name))}

  fits heterogeneous variances across \code{Location} for the random effect
  \code{Name}, while

  \code{random = ~ vsm(dsm(Environment), ar1m(Row), ar1m(Column),
  ism(Name), Gu=Ainv)}

  specifies an arbitrary Kronecker product of covariance structures before the
  relationship matrix for \code{Name}.

  Covariance constructors currently available include:

  \code{\link{ism}} for identity covariance;

  \code{\link{dsm}} and \code{\link{atm}} for diagonal and selected-level
  diagonal covariance structures;

  \code{\link{usm}} for an unstructured covariance matrix;

  \code{\link{csm}} and \code{\link{corgm}} for correlation structures;

  \code{\link{ar1m}}, \code{\link{ar2m}}, and \code{\link{ar3m}} for
  autoregressive covariance structures;

  \code{\link{mam}} (including \code{ma1m} and \code{ma2m}) for moving-average
  covariance structures;

  \code{\link{toeplitzm}} for a general positive-definite Toeplitz correlation
  structure;

  \code{\link{fam}} for factor-analytic covariance;

  \code{\link{rrcm}} for the reduced-rank covariance approximation;

  \code{\link{antem}} for antedependence covariance;

  \code{\link{maternm}} for Matern spatial covariance;

  \code{\link{sar}} and \code{\link{car}} for simultaneous and conditional
  autoregressive spatial covariance models; and

  \code{\link{ownm}} for a fixed or user-defined covariance function.

  Design-matrix utilities such as \code{\link[enhancer]{overlay}},
  \code{\link{spl2Dc}}, \code{\link[enhancer]{leg}}, and
  \code{\link[enhancer]{redmm}} can still be combined with \code{vsm()} when
  appropriate.

  See \code{\link{vsm}} and the individual covariance-constructor help pages
  for details on parameterizations and examples.

  A single relationship term may request Lee--van der Werf rotation with
  \code{vsm(..., Gu=Gu, rotation=TRUE)}. With \code{henderson=TRUE}, random
  coefficients are represented in the eigenbasis and the diagonal relationship
  precision is used in the mixed model equations. With \code{henderson=FALSE},
  the equivalent transformed observation-covariance parameterization is used.
  Public random effects, fitted values, and residuals are returned in the
  original basis. The current implementation requires a complete balanced
  Gaussian layout, identity residual covariance, default identity \code{W},
  one rotated relationship term, and \code{computeCi=0} during fitting.

  }
  \item{rcov}{A formula specifying the \strong{residual covariance structure}.
  The default is \emph{rcov = ~ units}, which fits a homogeneous residual
  variance.

  Structured residual covariance is specified with the same
  \code{\link{vsm}} and CovarianceFactor constructors used for random effects.
  The residual term should be expressed as a single \code{vsm()} term; arbitrary
  covariance complexity is represented by placing multiple covariance factors
  inside that term rather than by summing multiple residual terms.

  The final term is normally \code{ism(units)}, where \code{units} is the
  observation-level incidence generated internally by \code{mmes}. For example,

  \code{rcov = ~ vsm(dsm(Environment), ism(units))}

  fits heterogeneous residual variances across environments, while

  \code{rcov = ~ vsm(ar1m(Row), ar1m(Column), ism(units))}

  fits a separable row-by-column AR1 residual covariance. More generally,
  structures such as \code{\link{usm}}, \code{\link{csm}},
  \code{\link{toeplitzm}}, \code{\link{maternm}}, \code{\link{sar}},
  \code{\link{car}}, and \code{\link{ownm}} can be used when their incidence
  layout is appropriate for the residual model.

  Residual covariance factors must identify exactly one covariance-product
  coordinate for each observation. \code{mmes()} constructs the residual block
  and local-coordinate indexing internally, so the data do \strong{not} need to
  be manually sorted by the variables defining the residual structure.

  As with random effects, \code{vsm()} owns one residual variance scale
  \eqn{\sigma^2}; the preceding covariance constructors provide dimensionless
  covariance shapes that are combined through an arbitrary Kronecker product.

  See \code{\link{vsm}} and the individual covariance-constructor help pages
  for details.

  }
  \item{data}{Optional data frame containing variables used by the model.
  Model expressions are evaluated using normal R scoping rules: variables are first
  sought in \code{data} and may also be obtained from the environments associated
  with the model formulas or from the calling environment. Consequently, auxiliary
  objects such as relationship matrices, coordinate vectors, or grouping variables
  do not have to be copied into \code{data}. If \code{data} is omitted, variables
  are evaluated from the formula/calling environment.

  Observation-level variables used by the response, fixed effects, random effects,
  and residual covariance specification must nevertheless have compatible lengths.
  Internally, \code{mmes()} constructs a common observation mask and applies it
  consistently to all observation-level model components.}
  \item{W}{Weights matrix (e.g., when covariance among plots exist). Internally
  W is squared and inverted as Wsi = solve(chol(W)), then the residual matrix
  is calculated as R = Wsi*O*Wsi.t(), where * is the matrix product, and O 
  is the original residual matrix.}
   \item{nIters}{Maximum number of iterations allowed.}
  \item{tolParConvLL}{Convergence criteria based in the change of log-likelihood
  between iteration i and i-1.}
  \item{tolParConvNorm}{When using the Henderson method this argument is the 
  convergence criteria based in the norm proposed by Jensen, Madsen and
  Thompson (1997):

  e1 = || InfMatInv.diag()/sqrt(N) * dLu ||

  where InfMatInv.diag() is the diagonal of the inverse of the information
  matrix, N is the total number of variance components, and dLu is the vector 
  of first derivatives.

  }
  \item{tolParInv}{Tolerance parameter for matrix inverse used when singularities 
  are encountered in the estimation procedure. By default the value is 1e-06. This
  parameter should be fairly small because the it is used to bend matrices like the
  information matrix in the henderson algorithm or the coefficient matrix when it
  is not positive-definite.}
   \item{naMethodX}{Missing-data policy for variables entering the fixed-effect
   design matrix. The default, \code{"exclude"}, marks observations with missing
   fixed-effect information for removal from the common observation set.
   \code{"include"} retains the historical sommer behavior for supported cases by
   imputing missing fixed-effect covariates. This option should be used carefully,
   because imputation changes the fitted design matrix.}

   \item{naMethodY}{Missing-data policy for the response variable(s). The default,
   \code{"exclude"}, removes observations with missing response information from the
   common observation set. The historical \code{"include"} and \code{"include2"}
   behaviors are retained where supported: \code{"include"} imputes missing responses,
   whereas \code{"include2"} is intended for multivariate-response situations in which
   an observation may be retained when at least one response is observed.}

   \item{naMethodRandom}{Missing-data policy for observation-level variables used to
   construct random-effect terms. The default is \code{"exclude"}. Under this policy,
   an observation for which a required random-effect classification or covariate cannot
   be evaluated is excluded from the common observation set. Random terms are evaluated
   before the final observation subset is applied, so the same inclusion mask is used
   consistently for the response, fixed-effect design, and all random-effect design
   matrices.}

   \item{naMethodR}{Missing-data policy for observation-level variables used to define
   the residual covariance structure. The default is \code{"exclude"}. Observations
   lacking information required to assign their residual covariance-product coordinate
   are excluded from the common observation set. The residual block and local-coordinate
   layout are then constructed on the retained observations.}

   \item{returnParam}{A TRUE/FALSE value to indicate if the program should 
   return the parameters to be used for fitting the model instead of fitting
   the model.}
   \item{dateWarning}{A TRUE/FALSE value to indicate if the program should warn 
   you when is time to update the sommer package.}
  \item{verbose}{A TRUE/FALSE value to indicate if the program should return 
  the progress of the iterative algorithm.}
  \item{stepWeight}{A vector of values (of length equal to the number of 
  iterations) indicating the weight used to multiply the update (delta) for
  variance components at each iteration. If NULL the 1st iteration will be 
  multiplied by 0.5, the 2nd by 0.7, and the rest by 0.9. This argument can
  help to avoid that variance components go outside the parameter space in the 
  initial iterations which happens very often with the AI method but it can be 
  detected by looking at the behavior of the likelihood. In that case you may
  want to give a smaller weight.}
  \item{emWeight}{A vector of values (of length equal to the number of iterations)
  indicating with values between 0 and 1 the weight assigned to the EM
  information matrix; the values 1 - emWeight are applied to the AI information
  matrix to produce the joint information matrix used in the variance-component
  update. This creates an EM-heavy warm start early in the fit and gradually
  transitions to an AI-dominated update as the optimizer approaches convergence.
  The default schedule is a smooth decay from 1 at the first iteration to 0.05
  by the last iteration, i.e. the fit begins in a strongly stabilized EM regime
  and then moves toward the classical AI update. Values outside \code{[0, 1]} are
  rejected.}
  \item{contrasts}{an optional list. See the contrasts.arg of model.matrix.default.}
  \item{getPEV}{A logical value indicating whether prediction error variance
  results should be organized and returned when they are requested through
  \code{computeCi}.}
  \item{henderson}{A logical value indicating which REML engine to use. \code{TRUE}
  (the default) uses the Henderson mixed-model-equations algorithm, efficient when
  there are many more records than coefficients to estimate. \code{FALSE} uses a
  direct-inversion engine that inverts the phenotypic covariance matrix in
  observation space instead; it is less efficient than \code{henderson=TRUE} when
  there are many more records than coefficients, but can be more efficient when
  there are many more coefficients to estimate than records available (e.g.,
  marker/SNP-BLUP-style models). Both engines share the same \code{vsm()}
  covariance-structure parameterization. The direct-inversion engine currently
  supports \code{computeCi} values \code{0} and \code{2} only (no Takahashi
  selected-inverse mode), has no \code{solver}/PCG/CHOLMOD options, and requires a
  single response column (use the long-format \code{vsm(usm(trait), ...)}
  convention for multi-trait models, as with the Henderson engine).}
    \item{computeCi}{An integer indicating whether post-fit prediction error variance information
  and/or the inverse of the mixed model coefficient matrix should be computed when the
  Henderson algorithm is used (\code{henderson=TRUE}). The available options are:

  \code{0}: no additional inverse-related computation is performed after convergence.
  The full inverse of the coefficient matrix (\code{Ci}) is not formed and
  \code{uPevList} is not computed. This is the fastest option and is the default.

  \code{1}: prediction error variances are computed using the Takahashi sparse inverse
  subset algorithm. This method uses the sparse LDLT factorization of the coefficient
  matrix and recursively evaluates only those elements of the inverse that belong to
  the sparsity pattern induced by the factorization. In particular, the diagonal
  elements required for prediction error variances are obtained without constructing
  the complete inverse matrix. Therefore, \code{uPevList} is returned while \code{Ci}
  is not fully materialized.

  With \code{solver="cholmod"} (see below), the supernodal factorization used during
  REML iterations has no Takahashi selected-inverse equivalent. In that case
  \code{computeCi=1} automatically performs one additional sparse LDLT factorization
  of the converged coefficient matrix after REML convergence, purely to extract the
  selected inverse subset. This is a single one-time cost, not repeated every
  iteration, so the iterative speed benefit of \code{solver="cholmod"} is retained.

  \code{2}: the complete inverse of the coefficient matrix is computed. The prediction
  error variances in \code{uPevList} are then extracted from the appropriate diagonal
  elements of \code{Ci}. This option is the most computationally and memory intensive,
  especially for large mixed model equation systems. With \code{solver="cholmod"} the
  full inverse is obtained directly from the supernodal factor (no extra LDLT
  refactorization is needed for this mode).

  For large models, \code{computeCi=0} is recommended during model fitting. If prediction
  error variances or the full inverse are needed afterwards, they can be obtained with
  \code{\link{postPEV}} without refitting the model.}

   \item{solver}{Linear-system solver used by the Henderson implementation.
   \code{"auto"} (default) picks a solver automatically based on the density of the
   random-effect relationship matrices (\code{Gu}) actually supplied: if any random
   effect uses a dense \code{Gu} (e.g. a genomic or marker-based relationship matrix,
   which is typically close to fully dense), \code{"cholmod"} is selected; otherwise
   (e.g. sparse pedigree-based relationship matrices, or no \code{Gu} at all)
   \code{"ldlt"} is selected. Passing an explicit value (\code{"ldlt"}, \code{"pcg"},
   or \code{"cholmod"}) always overrides the automatic choice.
   \code{"ldlt"} uses the sparse direct simplicial LDLT factorization path.
   \code{"pcg"} selects the iterative preconditioned conjugate-gradient path.
   \code{"cholmod"} uses the supernodal (blocked, BLAS-3) sparse Cholesky
   factorization bundled with R's \code{Matrix} package (SuiteSparse CHOLMOD),
   accessed through \code{Matrix}'s public C API - no additional system library needs
   to be installed, since \code{Matrix} is already a dependency of \code{sommer}.
   Supernodal factorization can be faster than the simplicial LDLT path for large
   coefficient matrices (e.g., models built on large pedigree or genomic relationship
   matrices), because it reorganizes the sparse factorization into dense blocks that
   are handled with threaded BLAS-3 kernels rather than column-by-column updates; the
   speed-up depends on the BLAS linked to the running R installation. \code{solver="cholmod"}
   currently supports all \code{computeCi} values (see above for how \code{computeCi=1}
   is handled for this solver).
   The covariance-model specification is independent of this choice: all solvers
   receive the same generic CovarianceFactor descriptors produced by \code{vsm()}.}

   \item{pcgTol}{Convergence tolerance used by the PCG linear solver. The default is
  \code{1e-8}. This argument is used when \code{solver="pcg"}. Smaller values request
   more accurate iterative solves but may increase computation time.}

   \item{pcgMaxIters}{Maximum number of PCG iterations for an individual linear solve.
   A value of \code{0} lets the C++ solver choose its internal iteration limit from the
   problem dimension. This argument is used only by the PCG path.}

   \item{pcgTraceProbes}{Number of stochastic probe vectors used by the PCG-based
   machinery for trace quantities required during REML calculations. The default is
  \code{8}. Increasing this value can reduce stochastic approximation variability at
   additional computational cost. This argument is relevant to \code{solver="pcg"}.}

   \item{pcgLanczosSteps}{Number of Lanczos steps used by the PCG-based stochastic
  log-determinant calculations. The default is \code{20}. Larger values can improve
   the approximation at additional computational cost. This argument is relevant to
   \code{solver="pcg"}.}

   \item{REML}{Logical, \code{TRUE} (default) fits variance/covariance parameters by
   restricted maximum likelihood (REML), the classical choice for BLUP-style variance
   component estimation. Setting \code{REML=FALSE} fits by maximum likelihood (ML)
   instead: the log-likelihood, score equations, and reported \code{monitor}/log-lik
   values drop the REML "restriction" term that accounts for the degrees of freedom
   used to estimate the fixed effects. ML estimates of variance components are known
   to be biased downward relative to REML, but unlike REML log-likelihoods, ML
   log-likelihoods are directly comparable (via a likelihood-ratio test) across models
   that differ in their \strong{fixed}-effects specification - REML likelihoods are only
   comparable across models that share the same fixed effects. \code{REML=FALSE}
   currently requires \code{solver="ldlt"} or \code{solver="cholmod"} (\code{solver=
   "auto"} already resolves to one of these). \code{solver="pcg"} is not yet supported
   with \code{REML=FALSE}. Prediction error variances (\code{uPevList}/\code{computeCi})
   are still computed from the REML-style inverse regardless of \code{REML}.}

  \item{family}{A \code{stats::family()} object. The default Gaussian identity
  family uses the ordinary linear mixed-model implementation. Other families are
  fitted by penalized quasi-likelihood (PQL): each outer iteration forms an IRLS
  working response and weights, then fits the resulting weighted Gaussian mixed
  model with the Henderson solver. The initial implementation supports one numeric
  response. Binomial and Poisson working residual dispersions are fixed to one.
  Formula offsets are included on the link scale. A supplied symmetric positive
  definite \code{W} is combined with the IRLS precision after observation filtering,
  and may be sparse and non-diagonal.}

  \item{pqlControl}{A named list controlling non-Gaussian PQL fits. Supported
  entries are \code{maxit} (outer-iteration limit, default 20) and \code{tol}
  (relative deviance convergence tolerance, default \code{1e-5}).}

  \item{.pqlInner}{Internal-use logical indicating that \code{mmes()} is being
  called for a Gaussian working-model fit within the PQL iteration. Users
  should not set this argument directly.}

  \item{.pqlFixedDispersion}{Internal-use logical indicating that the working
  residual dispersion is fixed, as required for binomial and Poisson PQL
  models. Users should not set this argument directly.}

  \item{.pqlWorkingPrecision}{Internal-use vector containing the current IRLS
  working precisions after observation filtering. Users should not set this
  argument directly.}

  \item{.pqlBaseW}{Internal-use base observation precision matrix supplied by
  the outer PQL fit before combination with the IRLS working precisions. Users
  should not set this argument directly.}

  \item{.pqlBaseFactor}{Internal-use cached factor of \code{.pqlBaseW}, reused
  across PQL iterations. Users should not set this argument directly.}
}

\details{

\strong{Evaluation environments and observation filtering}

The current interface follows R-style model evaluation. Variables referenced by
\code{fixed}, \code{random}, and \code{rcov} may be supplied as columns of
\code{data} or resolved from the formula/calling environment. This makes it possible,
for example, to keep relationship matrices, external grouping variables, or spatial
coordinates outside the analysis data frame.

Missingness is handled through a common observation-selection mechanism. The response,
fixed effects, random effects, and residual covariance specification are evaluated
against the original observation set and their validity is combined into a single
inclusion mask. The same retained rows are then used for the response vector,
fixed-effect design matrix, random-effect design matrices, optional weight matrix,
and residual covariance indexing. This avoids independently subsetting different
parts of the mixed model.

The fitted object contains \code{obsInfo}, an observation-level audit table describing
the filtering process. It records the original row index, validity of the response,
fixed, random, and residual components, whether the observation was retained, and,
when applicable, the reason for exclusion. This can be useful for diagnosing model
specifications involving missing values.

\strong{Estimated marginal means}

When the suggested \pkg{emmeans} package is installed, univariate \code{mmes}
fits can be supplied directly to \code{emmeans::emmeans()}. The resulting means
are population-level fixed-effect marginal means with covariance obtained from
the mixed-model equations. Inference uses asymptotic degrees of freedom
(\code{df = Inf}). Multivariate-response fits are not currently supported by
this interface.

Random and residual formulas are evaluated as R language objects rather than by
splitting their textual representation at plus signs. Consequently, arithmetic or
nested expressions inside model terms are not interpreted as separate top-level
random or residual terms.

The use of this function requires a good understanding of mixed models. 
Please review the 'sommer.quick.start' vignette and pay attention to details 
like format of your random and fixed variables (e.g. character and factor 
variables have different properties when returning BLUEs or BLUPs).

\strong{For tutorials} on how to perform different analysis with sommer please
look at the vignettes by typing in the terminal:

vignette("v1.sommer.quick.start")

vignette("v2.sommer.changes.and.faqs")

vignette("v3.sommer.qg")

vignette("v4.sommer.gxe")

\strong{Citation}

Type \emph{citation("sommer")} to know how to cite the sommer package in your publications.

\strong{Special variance structures}

 \code{\link{vsm}(\link{atm}(x,levels),ism(y))}

 can be used to specify heterogeneous variance for the "y" covariate at specific
 levels of the covariate "x", e.g., \emph{random=~vsm(at(Location,c("A","B")),ism(ID))}
 fits a variance component for ID at levels A and B of the covariate Location.

 \code{\link{vsm}(\link{dsm}(x),ism(y))}

 can be used to specify a diagonal covariance structure for the "y" covariate 
 for all levels of the covariate "x", e.g., \emph{random=~vsm(dsm(Location),ism(ID))} 
 fits a variance component for ID at all levels of the covariate Location.

 \code{\link{vsm}(\link{usm}(x),ism(y))}

 can be used to specify an unstructured covariance structure for the "y" covariate 
 for all levels of the  covariate "x", e.g., \emph{random=~vsm(usm(Location),ism(ID))}  
 fits variance and covariance components for ID at all levels of the covariate Location.
 
 \code{\link{vsm}(\link{usm}(\link[enhancer]{rrm}(x,y,z,nPC)),ism(y))}

 can be used to specify an unstructured covariance structure for the "y" effect
 for all levels of the  covariate "x", and a response variable "z", e.g., 
 \emph{random=~vsm(rrm(Location,ID,response, nPC=2),ism(ID))}  fits a reduced-rank 
 factor analytic covariance for ID at 2 principal components of the covariate Location.

  \code{\link{vsm}(ism(\link[enhancer]{overlay}(...,rlist=NULL,prefix=NULL)))}

 can be used to specify overlay of design matrices between consecutive random 
 effects specified, e.g., \emph{random=~vsm(ism(overlay(male,female)))} overlays
 (overlaps) the incidence matrices for the male and female random effects to obtain
 a single variance component for both effects. The `rlist` argument is a list with 
 each element being a numeric value that multiplies the incidence matrix to be 
 overlayed. See \code{\link[enhancer]{overlay}} for details.Can be combined with vsm().
 
   \code{\link{vsm}(ism(\link[enhancer]{redmm}(x,M,nPC)))}

 can be used to create a reduced model matrix of an effect (x) assumed to be a 
 linear function of some feature matrix (M), e.g., \emph{random=~vsm(ism(redmm(x,M)))} 
 creates an incidence matrix from a very large set of features (M) that belong 
 to the levels of x to create a reduced model matrix. See \code{\link[enhancer]{redmm}} 
 for details.Can be combined with vsm().

  \code{\link{vsm}(\link[enhancer]{leg}(x,n),ism(y))}

  can be used to fit a random regression model using a numerical variable \code{x} 
  that marks the trayectory for the random effect \code{y}. The leg function can 
  be combined with the special functions \code{dsm}, \code{usm} \code{at} and \code{csm}. 
  For example \emph{random=~vsm(leg(x,1),ism(y))} or \emph{random=~vsm(usm(leg(x,1)),ism(y))}.

  \code{\link{spl2Dc}(x.coord, y.coord, at.var, at.levels))}

  can be used to fit a 2-dimensional spline (e.g., spatial modeling) using 
  coordinates \code{x.coord} and \code{y.coord} (in numeric class) assuming 
  multiple variance components. The 2D spline can be fitted at specific 
  levels using the \code{at.var} and \code{at.levels} arguments. For example
  \emph{random=~spl2Dc(x.coord=Row,y.coord=Range,at.var=FIELD)}.

\strong{Covariance between random effects}

 \code{\link{covm}( \link{vsm}(\link{ism}(ran1)), \link{vsm}(\link{ism}(ran2)) )}

 can be used to specify covariance between two different random effects, e.g.,
 \emph{random=~covm( vsm(ism(x1)), vsm(ism(x2)) )} where two random effects in 
 their own vsm() structure are encapsulated. Only applies for simple random effects.
 
\strong{S3 methods}

S3 methods are available for some parameter extraction such as
\code{\link{fitted.mmes}}, \code{\link{residuals.mmes}}, \code{\link{summary.mmes}},
\code{\link{randef}}, \code{\link{coef.mmes}}, \code{\link{anova.mmes}},
\code{\link{plot.mmes}}, and \code{\link{predict.mmes}} to obtain adjusted means.
In addition, the \code{\link{vpredict}} function (replacement of the pin 
function) can be used to estimate standard errors for linear combinations of 
variance components (e.g., ratios like h2). The \code{\link{r2}} function calculates
reliability. 


\strong{Additional Functions}

Additional functions for genetic analysis have been included such as relationship 
matrix building (\code{\link{A.mat}}, \code{\link{D.mat}}, \code{\link{E.mat}},
\code{\link{H.mat}}),  build a genotypic hybrid marker matrix (\code{\link[enhancer]{build.HMM}}), 
plot of genetic maps (\code{\link[enhancer]{map.plot}}), and manhattan plots
(\code{\link[enhancer]{manhattan}}). If you need to build a pedigree-based relationship 
matrix use the \code{getA} function from the pedigreemm package.

\strong{Bug report and contact}

If you have any technical questions or suggestions please post it in
https://stackoverflow.com or https://stats.stackexchange.com

If you have any bug report please go to https://github.com/covaruber/sommer or 
send me an email to address it asap, just make sure you have read the vignettes 
carefully before sending your question.

\strong{Example Datasets}

The package has been equiped with several datasets to learn how to use the sommer package:

* \code{\link[enhancer]{DT_halfdiallel}}, \code{\link[enhancer]{DT_fulldiallel}} and \code{\link[enhancer]{DT_mohring}} 
datasets have examples to fit half and full diallel designs.

* \code{\link[enhancer]{DT_h2}} to calculate heritability

* \code{\link[enhancer]{DT_cornhybrids}} and \code{\link[enhancer]{DT_technow}} datasets to perform
genomic prediction in hybrid single crosses

* \code{\link[enhancer]{DT_wheat}} dataset to do genomic prediction in single crosses in 
species displaying only additive effects.

* \code{\link[enhancer]{DT_cpdata}} dataset to fit genomic prediction models within a
biparental population coming from 2 highly heterozygous parents including additive, 
dominance and epistatic effects.

* \code{\link[enhancer]{DT_polyploid}} to fit genomic prediction and GWAS analysis in polyploids.

* \code{\link[enhancer]{DT_gryphon}} data contains an example of an animal model including 
pedigree information.

* \code{\link[enhancer]{DT_btdata}} dataset contains an animal (birds) model.

* \code{\link[enhancer]{DT_legendre}} simulated dataset for random regression model.

* \code{\link[enhancer]{DT_sleepstudy}} dataset to know how to translate lme4 models to 
sommer models.

* \code{\link[enhancer]{DT_ige}} dataset to show how to fit indirect genetic effect models.

\strong{Models Enabled}

For details about the models enabled and more information about the covariance 
structures please check the help page of the package (\code{\link{sommer}}).

}

\value{
If all parameters are correctly indicated the program will return a list with the following information:

\item{data}{the dataset used in model fitting after application of the common
              observation mask.}
\item{obsInfo}{an observation-level audit table describing the common missing-data
               filter. It includes the original row index, validity indicators for
               response, fixed, random, and residual model components, the final
               inclusion indicator, and the reason for exclusion when applicable.}
\item{Dtable}{the table to be used for the predict function to help the program
              recognize the factors available.}
\item{llik}{the vector of log-likelihoods across iterations}
\item{b}{the vector of fixed effect.}
\item{u}{the vector of random effect.}
\item{bu}{the vector of fixed and random effects together.}
\item{Ci}{the inverse of the coefficient matrix.}
\item{Ci_11}{the inverse of the coefficient matrix pertaining to the fixed effects.}
\item{theta}{a list of estimated variance covariance matrices. Each element of 
             the list corresponds to the different random and residual components}
\item{covPar}{the fitted covariance parameters in their descriptor-defined
reported coordinates.}
\item{covParNative}{a data frame containing all fixed and estimated covariance
parameters in the native model-facing scale defined by each CovarianceFactor.
For example, \code{dsm()} entries are level-specific variances rather than
variance ratios. This is the same result returned by \code{covparams_mmes()}.}
\item{covParNativeSE}{the native-scale covariance-parameter table augmented
with delta-method standard errors and Z ratios. This is the same result
returned by \code{covparams_mmes_se()}.}
\item{theta_se}{inverse of the information matrix.}
\item{InfMat}{information matrix.}
\item{monitor}{The values of the variance-covariance components across iterations 
              during the REML estimation.}
\item{AIC}{Akaike information criterion}
\item{BIC}{Bayesian information criterion}
\item{convergence}{a TRUE/FALSE statement indicating if the model converged.}
\item{partitions}{a list where each element contains a matrix indicating where 
                  each random effect starts and ends.}
\item{partitionsX}{a list where each element contains a matrix indicating where 
                   each fixed effect starts and ends.}
\item{percDelta}{the matrix of percentage change in deltas (see tolParConvNorm argument).}
\item{normMonitor}{the matrix of the three norms calculated (see tolParConvNorm argument).}
\item{toBoundary}{the matrix of variance components that were forced to the
                   boundary across iterations.}
\item{Cchol}{the Cholesky decomposition of the coefficient matrix.}
\item{y}{the response vector.}
\item{W}{the column binded matrix W = [X Z] }
\item{uList}{a list containing the BLUPs in data frame format where rows are 
             levels of the random effects and column the different factors at 
             which the random effect is fitted. This is specially useful for 
             diagonal and unstructured models.}
\item{uPevList}{a list containing the BLUPs in data frame format where rows are 
                levels of the random effects and column the different factors at 
                which the random effect is fitted. This is specially useful for 
                diagonal and unstructured models.}
\item{args}{the fixed, random and residual formulas from the mmes model.}
\item{constraints}{The vector of constraints.}

}
\references{

Covarrubias-Pazaran G. Genome assisted prediction of quantitative traits using 
the R package sommer. PLoS ONE 2016, 11(6): doi:10.1371/journal.pone.0156744

Jensen, J., Mantysaari, E. A., Madsen, P., and Thompson, R. (1997). Residual 
maximum likelihood estimation of (co) variance components in multivariate mixed 
linear models using average information. Journal of the Indian Society of
Agricultural Statistics, 49, 215-236.

Sanderson, C., & Curtin, R. (2025). Armadillo: An Efficient Framework for 
Numerical Linear Algebra. arXiv preprint arXiv:2502.03000.

Gilmour et al. 1995. Average Information REML: An efficient algorithm for 
variance parameter estimation in linear mixed models. Biometrics 51(4):1440-1450.

}
\author{

Coded by Giovanny Covarrubias-Pazaran with contributions of Christelle Fernandez 
Camacho, Johan Aparicio-Arce, and Claudio Flavio 5.1 jaja (laughing in Spanish).

}
\examples{

####=========================================####
#### For CRAN time limitations most lines in the
#### examples are silenced with one '#' mark,
#### remove them and run the examples
####=========================================####

data(DT_example, package="enhancer")
DT <- DT_example
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

\donttest{

####===========================================####
#### Univariate heterogeneous variance models  ####
####===========================================####
DT=DT[with(DT, order(Env)), ]
## Compound simmetry (CS) + Diagonal (DIAG) model
ans2 <- mmes(Yield~Env,
             random= ~Name + vsm(dsm(Env),ism(Name)),
             rcov= ~ vsm(dsm(Env),ism(units)),
             data=DT)
summary(ans2)

####===========================================####
####  Univariate unstructured variance models  ####
####===========================================####

ans3 <- mmes(Yield~Env,
             random=~ vsm(usm(Env),ism(Name)),
             rcov=~vsm(dsm(Env),ism(units)),
             data=DT)
summary(ans3)


}

}
