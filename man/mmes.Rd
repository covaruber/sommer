\name{mmes}
\alias{mmes}
\title{

\strong{m}ixed \strong{m}odel \strong{e}quations \strong{s}olver

}
\description{

Fits linear mixed models by restricted maximum likelihood (REML) or maximum
likelihood (ML), and generalized linear mixed models by penalized
quasi-likelihood (PQL). Random and residual covariance structures are specified
with \code{\link{vsm}} and covariance constructors. Gaussian fits can use
Henderson mixed-model equations or direct inversion in observation space;
\code{solver} selects LDLT, CHOLMOD or PCG only for the Henderson engine.
Alternatively, \code{solveOnly=TRUE} solves the mixed-model equations at known
covariance parameters without estimating them. Numerical engines are implemented
in C++ using Armadillo and Eigen.


}

\usage{

mmes(fixed, random, rcov, data, W, weights=NULL,
     nIters=30, tolParConvLL=1e-04,
     tolParConvNorm=1e-04, tolParInv=1e-06,
     naMethodX="exclude", naMethodY="exclude",
     naMethodRandom="exclude", naMethodR="exclude",
     returnParam=FALSE, dateWarning=TRUE,
     verbose=TRUE, stepWeight=NULL, emWeight=NULL,
     contrasts=NULL, getPEV=TRUE, henderson=TRUE,
    computeCi=0, solver="auto", pcgTol=1.0e-8,
    pcgMaxIters=0, pcgTraceProbes=8,
    pcgLanczosSteps=20, REML=TRUE, vcc=NULL,
    family=stats::gaussian(), pqlControl=list(),
    .pqlInner=FALSE, .pqlFixedDispersion=FALSE,
    .pqlWorkingPrecision=NULL, .pqlBaseW=NULL,
    .pqlBaseFactor=NULL, acceleration="none", .pqlStart=NULL,
     factorScoreAugmentation="none",
     .factorScoreParameters=NULL,
    pcgPreconditioner="diagonal", pcgNystromRank=32L,
    solveOnly=FALSE, covPar=NULL)
}
\arguments{

  \item{fixed}{A formula specifying the \strong{response variable(s)} \strong{and fixed effects}, i.e:

  \emph{response ~ covariate}. Use one numeric response column; multi-trait
  models should be stacked into long format with \code{\link{stackTraits}},
  with the trait factor included in the fixed, random and residual formulas
  as appropriate. Formula offsets are supported.

  }

  \item{random}{A formula specifying the \strong{random effects}, e.g.,
  \emph{random = ~ genotype + year}. Omit this argument for a model without
  random effects. A simple random term that is not wrapped
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

  \code{\link{rrm}} for the reduced-rank covariance approximation;

  \code{\link{antem}} for antedependence covariance;

  \code{\link{maternm}} for Matern spatial covariance;

  \code{\link{metricm}} for exponential, power, gaussian, spherical and
  circular distance-based correlations (isotropic or anisotropic);

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
  Gaussian layout, one rotated relationship term, and \code{computeCi=0} during
  fitting. Residual structures are allowed when they are rotation invariant:
  residuals may be structured across rotation blocks (e.g.,
  \code{rcov=~vsm(usm(trait), ism(units))} or
  \code{rcov=~vsm(dsm(Env), ism(units))}), but must be independent and
  identically structured across the levels of the rotated relationship term.
  Residual structures that correlate observations of different levels (e.g.,
  \code{ar1m()}, \code{maternm()}, \code{sar()}, or \code{car()} over plots),
  residual variances that differ among levels within a rotation block, and
  \code{W} matrices that are not diagonal with a constant weight within each
  rotation block are rejected, as are non-Gaussian PQL fits.

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

  fits a separable row-by-column AR1 residual covariance, with one correlation
  per direction shared by every other Kronecker factor. For section-specific
  variances and correlations, use

  \code{rcov = ~ dsumm(vsm(ar1m(Row), ar1m(Column), ism(units)), by = Trial)}

  (see \code{\link{dsumm}}). More generally,
  structures such as \code{\link{usm}}, \code{\link{csm}},
  \code{\link{toeplitzm}}, \code{\link{maternm}}, \code{\link{sar}},
  \code{\link{car}}, and \code{\link{ownm}} can be used when their incidence
  layout is appropriate for the residual model.

  Residual covariance factors must identify exactly one covariance-product
  coordinate for each observation. \code{mmes()} constructs the residual block
  and local-coordinate indexing internally, so the data do \strong{not} need to
  be manually sorted by the variables defining the residual structure.

  \strong{Residual pairing key (multi-trait models).} For a long-format
  multi-trait model the residual must know which records belong to the same
  experimental unit. With \code{ism(units)} records are paired implicitly by
  their order of appearance within each trait; a warning is issued when
  repeated residual coordinates must be paired by row order. An explicit key is given by replacing
  \code{units} with a column identifying the unit, e.g.

  \code{rcov = ~ vsm(usm(trait), ism(record))}

  where \code{record} is the column created by \code{\link{stackTraits}}.
  Records sharing a key value form one residual block (missing traits are
  simply absent from it); a trait may appear only once per key value.

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
  \item{W}{Weights matrix (e.g., when covariance among plots exists).
  Omit this argument for unweighted Gaussian fitting; it has no explicit
  \code{NULL} default, so omission is preferable to passing \code{W=NULL}.
  \code{W} must be symmetric positive definite, with one row and column per
  original observation; it is subset using the common observation mask.
  Internally Wsi = solve(chol(W)), then the residual matrix
  is calculated as R = Wsi*O*Wsi.t(), where * is the matrix product, and O 
  is the original residual matrix. Weights can be combined with any residual
  structure, including Kronecker (multi-trait) residuals.}
  \item{weights}{Optional one-sided formula describing independent row blocks
  of the supplied \code{W} matrix, for example \code{weights = ~trial}. Rows
  with different formula groups must have zero cross-group entries in \code{W}.
  This is structural metadata only: it does not generate numeric weights.
  Declared blocks are processed independently by the Henderson weighting
  path; without this formula, \code{W} is handled without a known grouping.}
  \item{nIters}{Maximum number of covariance-estimation iterations, default
  30. In PQL this is the limit for each inner Gaussian fit; use
  \code{pqlControl$maxit} for the outer limit. It is not used by
  \code{solveOnly=TRUE}, which uses \code{pcgMaxIters}.}
  \item{tolParConvLL}{Absolute log-likelihood change tolerance between
  successive iterations, default \code{1e-4}. The Henderson engine stops
  when either this tolerance or \code{tolParConvNorm} is satisfied.}
  \item{tolParConvNorm}{When using the Henderson method this argument is the 
  convergence tolerance (default \code{1e-4}) based on the norm proposed by Jensen, Madsen and
  Thompson (1997):

  e1 = || InfMatInv.diag()/sqrt(N) * dLu ||

  where InfMatInv.diag() is the diagonal of the inverse of the information
  matrix, N is the total number of variance components, and dLu is the vector 
  of first derivatives. This criterion is not used by the direct engine.

  }
  \item{tolParInv}{Tolerance parameter for matrix inverse used when singularities 
  are encountered in the estimation procedure. By default the value is 1e-06. This
  parameter should be fairly small because it is used to bend matrices like the
  information matrix in the henderson algorithm or the coefficient matrix when it
  is not positive-definite.}
   \item{naMethodX}{Missing-data policy for variables entering the fixed-effect
   design matrix. The default, \code{"exclude"}, marks observations with missing
   fixed-effect information for removal from the common observation set.
  \code{"omit"} is an alias for \code{"exclude"}. \code{"fail"} stops
  if any required fixed-effect variable is missing. \code{"include"}
  (alias \code{"pass"}) retains rows at the filtering stage but does not
  impute values; downstream design construction may reject missing covariates.
  Impute explicitly before fitting if needed.}

   \item{naMethodY}{Missing-data policy for the response variable(s). The default,
   \code{"exclude"}, removes observations with missing response information from the
  common observation set. \code{"omit"} is an alias. \code{"fail"}
  stops on missing responses. \code{"include"} (alias \code{"pass"})
  bypasses filtering but does not impute responses, so missing values may
  prevent fitting. \code{"include2"} is not supported. For multi-trait
  models, stack the data and exclude only the rows with missing responses.}

   \item{naMethodRandom}{Missing-data policy for observation-level variables used to
   construct random-effect terms. The default is \code{"exclude"}. Under this policy,
   an observation for which a required random-effect classification or covariate cannot
   be evaluated is excluded from the common observation set. Random terms are evaluated
   before the final observation subset is applied, so the same inclusion mask is used
   consistently for the response, fixed-effect design, and all random-effect design
  matrices. \code{"omit"} is an alias for \code{"exclude"};
  \code{"fail"} stops on missing required values. \code{"include"}
  (alias \code{"pass"}) bypasses this filter without imputing classifications
  or covariates; the random-effect constructor must still produce a valid design.}

   \item{naMethodR}{Missing-data policy for observation-level variables used to define
   the residual covariance structure. The default is \code{"exclude"}. Observations
   lacking information required to assign their residual covariance-product coordinate
   are excluded from the common observation set. The residual block and local-coordinate
  layout are then constructed on the retained observations. \code{"omit"}
  is an alias for \code{"exclude"}; \code{"fail"} stops on missing
  required values. \code{"include"} (alias \code{"pass"}) bypasses this
  filter without imputing residual coordinates; valid covariance coordinates
  are still required. Section-specific missing-coordinate behavior is
  described in \code{\link{dsumm}}.}

  \item{returnParam}{Logical, default \code{FALSE}. For a Gaussian
  covariance-estimation call, \code{TRUE} returns prepared design matrices,
  covariance descriptors, the parameter-identification table \code{vcParams},
  and fitting controls instead of fitting. It is not a request for fitted
  covariance estimates. PQL and \code{solveOnly=TRUE} are dispatched before
  this option and do not return this setup list.}
  \item{dateWarning}{Logical, default \code{TRUE}. Print an update reminder
  when the installed package date is more than 90 days old.}
  \item{verbose}{Logical, default \code{TRUE}. Print iterative progress.
  Some setup messages are printed independently of this setting.}
  \item{stepWeight}{A vector of values (of length equal to the number of 
  iterations) indicating the weight used to multiply the update (delta) for
  covariance parameters at each iteration. If \code{NULL}, weights are 0.9,
  except that the first two iterations whose \code{emWeight} is at most 0.5
  use 0.5 and 0.7 when both exist. If fewer than two such iterations exist,
  the first two iterations use those weights instead. A supplied scalar is
  repeated; other supplied vectors are recycled to \code{nIters}.
  All weights must be finite and positive. This argument can
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
  The default schedule is a logarithmic decay from 1 to 0.03 over the first
  \code{min(nIters, 13)} iterations, followed by a 0.03 floor. For a
  one-iteration fit the default is 1. Supply one value per iteration when
  specifying a custom schedule.
  Values outside \code{[0, 1]} are
  rejected.}
  \item{contrasts}{Optional named list of contrasts for fixed-effect factors,
  passed as \code{contrasts.arg} to \code{\link[stats]{model.matrix}}.
  The default \code{NULL} uses the current R contrast options.}
  \item{getPEV}{A logical value indicating whether prediction error variance
  results should be organized and returned when they are requested through
  \code{computeCi}.}
  \item{henderson}{A logical value indicating which Gaussian REML/ML engine to use. \code{TRUE}
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
  The automatic threshold is density strictly greater than 0.2 in any
  prepared random-effect relationship precision.
   \code{"ldlt"} uses the sparse direct simplicial LDLT factorization path.
   \code{"pcg"} selects the iterative preconditioned conjugate-gradient path.
  For REML fits with \code{computeCi=0}, random-effect precisions produced by
  \code{vsm()} are applied matrix-free using their Kronecker factors, avoiding
  explicit assembly of \eqn{K^{-1} \otimes A^{-1}} during optimization. The
   complete coefficient matrix is materialized once after convergence for
   compatibility with prediction and post-fit methods. Other PCG configurations
   transparently use the assembled sparse path.
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
  Eligible coefficient systems use an exact dense block-Schur engine. Admission
  uses estimated dense storage rather than a fixed effect-count limit.
  \code{options(sommer.mme.denseMemoryMB=4096)} sets a 4096 decimal MB budget;
  the default is 1000 MB. This is not a bound on total process memory: sparse
  matrices, temporary workspaces and fitted output require additional memory.
  The budget also controls coefficient inverse-cache admission, with a separate
  workspace estimate. It does not change the residual engine's layout rules.
  For a coupled Kronecker random term, the full coupled dimension is used, not
  only the main-effect dimension. The fitted object's \code{engineDiagnostics}
  reports \code{denseMemoryMB} and \code{blockSchurEstimatedDenseMB}.
  For REML with one coupled random term whose coefficient-block graph consists
  of paths, dense local blocks and a sufficiently small border, an exact
  bordered block-chain engine is selected first. This can exploit AR1 random
  precision using main-effect-sized blocks without factoring the full coupled
  dense matrix. Unsupported graphs retain the usual dense/sparse dispatch.
  \code{engineDiagnostics$blockChainActive} identifies the chain path;
  \code{blockSchurActive} is true for either bordered engine. Residual layout
  heuristics are unchanged, and \code{computeCi=1} still performs its final
  one-time LDLT inverse-subset extraction.
  Eligible FA/RR random terms use an exact latent-factor Schur representation
  automatically, including free loadings and specific variances. This requires
  Gaussian Henderson REML, one covariance shaping factor in the selected term,
  well-conditioned positive specifics, dense relationship precision, disjoint
  covariance-coordinate incidence blocks and diagonal effective residual
  precision. Virtual factor equations are eliminated internally; the original
  coefficient ordering, marginal C, analytic AI updates and standard errors
  are preserved. No \code{factorScoreAugmentation} setting is needed.
  \code{engineDiagnostics$factorSchurActive} identifies this path and
  \code{factorSchurFallbacks} counts numerical latent-factor failures that
  reverted to the original marginal CHOLMOD solve. The dense-storage estimate
  includes virtual border equations, while total process memory also includes
  the explicitly assembled marginal C. Incompatible models keep existing paths.
   The covariance-model specification is independent of this choice: all solvers
   receive the same generic CovarianceFactor descriptors produced by \code{vsm()}.}

   \item{pcgTol}{Convergence tolerance used by the PCG linear solver. The default is
  \code{1e-8}. This argument is used when \code{solver="pcg"}. Smaller values request
   more accurate iterative solves but may increase computation time.}

   \item{pcgMaxIters}{Maximum number of PCG iterations for an individual linear solve.
   A value of \code{0} lets the C++ solver choose its internal iteration limit from the
  problem dimension. This argument is used by \code{solver="pcg"} and
  by the separate known-parameter solver when \code{solveOnly=TRUE}.}

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

  \item{vcc}{Optional data frame of equality/scaling constraints between
  covariance parameters (columns \code{parameter}, \code{group} and optionally
  \code{scale}), identified through the \code{vcParams} table returned by
  \code{returnParam=TRUE} and stored in the fit. See \code{\link{vcc}}.}

  \item{family}{A \code{stats::family()} object. The default Gaussian identity
  family uses the ordinary linear mixed-model implementation. A
  \code{\link{familym}} object specifies trait-specific families. Other families are
  fitted by penalized quasi-likelihood (PQL): each outer iteration forms an IRLS
  working response and weights, then fits the resulting weighted Gaussian mixed
  model with the Henderson solver. The initial implementation supports one numeric
  response. Any family with a valid link can be used, e.g. \code{binomial()}
  (logit, probit, cloglog), \code{poisson()}, \code{Gamma()},
  \code{inverse.gaussian()}, \code{gaussian(link="log")}, \code{quasipoisson()}
  and \code{MASS::negative.binomial(theta)}. Binomial counts out of a known
  number of trials are supplied as a single proportion response built with
  \code{\link{binm}}, whose trials are used as prior weights. Binomial,
  Poisson and negative binomial working residual dispersions are fixed to one;
  for other families the reported \code{dispersion} is the Pearson estimate.
  The negative binomial \eqn{\theta} is estimated (conditional on the random
  effects) with \code{pqlControl=list(estimateTheta=TRUE)} and returned in
  \code{theta.nb} and \code{theta.nb.se}. Because the working-model
  log-likelihood is not a GLMM likelihood, \code{llik}, \code{AIC} and
  \code{BIC} are \code{NA} for PQL fits and \code{anova()} likelihood ratio
  tests are refused; use \code{deviance} and Wald tests instead.
  Formula offsets are included on the link scale. A supplied symmetric positive
  definite \code{W} is combined with the IRLS precision after observation filtering,
  and may be sparse and non-diagonal.

  Long-format multi-trait data with a different family per trait are fitted by
  supplying \code{family = \link{familym}("trait", Disease = binomial(),
  Yield = gaussian())}. See \code{\link{familym}} for the required residual
  structure; in that case \code{dispersion} is a named vector with one value
  per trait.}

  \item{pqlControl}{A named list controlling non-Gaussian PQL fits. Supported
  entries are \code{maxit} (outer-iteration limit, default 20), \code{tol}
  (relative deviance convergence tolerance, default \code{1e-5}) and
  \code{estimateTheta} (estimate the negative binomial \eqn{\theta}, default
  \code{FALSE}) and \code{warmStart} (reuse covariance parameters from the
  preceding working fit when the observation and model layouts agree, default
  \code{TRUE}). The \code{pqlMonitor} also reports \code{warmStarted} and
  \code{innerIterations} and \code{symbolicAnalyses}. With the LDLT solver,
  compatible warm fits reuse symbolic factorizations and selected-inverse
  topology across PQL iterations while recomputing numeric factors. Other
  solvers reuse their symbolic factorizations within each working fit.}

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
  
  \item{acceleration}{Optional optimizer acceleration: \code{"none"} (default) or
  \code{"aitken"}. Aitken proposals are capped, applied every third accepted
  iteration during EM blending, and checked by the existing likelihood line
  search. Requires a deterministic Henderson solver; not available for PCG or
  direct inversion. The Henderson result's \code{engineDiagnostics} reports
  coefficient evaluations, line-search halvings, accelerated proposals,
  symbolic-analysis and ML selected-inverse topology counters, and active
  block-Schur partition dimensions, grouped weight-block counts and assembled
  PCG batch solve/iteration counts, plus peak sparse fill in the irregular
  residual fallback.}
  \item{.pqlStart}{Internal-use covariance warm-start state for successive
  PQL working fits. The default is \code{NULL}; users should not set it
  directly.}
  \item{factorScoreAugmentation}{FA/RR latent factor parameterization. The
  default \code{"none"} uses the marginal covariance MME. \code{"fixed-shape"}
  uses the augmented MME with fixed loading and specific-variance parameters;
  it is currently limited to Henderson REML, \code{computeCi=0}, and one FA/RR
  shaping factor per eligible random term. It replaces each \eqn{q d} marginal
  coefficient block with \eqn{(k+q)d} factor and specific-effect coefficients
  constrained to share one variance scale. \code{"profile"} estimates free
  FA/RR shape parameters with a bounded L-BFGS-B outer profile optimizer. It
  starts from the values in the constructors (so it does not first factor the
  expensive marginal FA/RR system); at each candidate shape it rebuilds the
  augmented design and refits shared random and residual scales by REML.
  The outer iteration budget is dimension-adjusted, with at least five and at
  most twelve outer iterations; \code{nIters} controls each inner AI fit.
  Inner covariance starts and symbolic factorizations are reused when the
  augmented layout matches; numeric factors are recomputed. The outer objective
  uses finite-difference gradients, is experimental, can require many augmented
  REML fits, and can converge to a local optimum; inspect
  \code{factorScoreProfile}. Shape-parameter standard errors are returned as
  \code{NA}; the inner AI information matrix is not their uncertainty estimate.

  Both augmented modes currently require Gaussian identity Henderson REML,
  \code{computeCi=0}, one FA/RR shaping factor per eligible random term, no
  rotation, and no user \code{vcc}. With \code{solver="auto"}, dense
  \code{Gu} selects CHOLMOD; its exact block-Schur path eliminates independent
  specific-effect blocks while retaining factor effects in the border. The
  returned \code{uList}, \code{bu}, fitted values, and prediction contrasts use
  the original effects; \code{C} remains augmented and is tagged with
  \code{CRepresentation="factor-score-augmented"}.}
  
  \item{.factorScoreParameters}{Internal-use list of FA/RR shape parameters
  on the working scale, indexed by random term. The profile optimizer uses
  these values to fix the shape parameters in each inner Gaussian fit.
  Users should not set this argument directly.}
  
  \item{pcgPreconditioner}{PCG preconditioner: \code{"diagonal"} (default) or
  \code{"nystrom"}. The Nystrom option builds a deterministic landmark
  approximation to assembled \eqn{C} and applies its Woodbury inverse as an SPD
  preconditioner; PCG continues to multiply by the exact \eqn{C} and verifies
  the true residual. Nystrom currently uses assembled-C PCG, so it disables the
  matrix-free shortcut.}
  
  \item{pcgNystromRank}{Positive integer number of deterministic landmarks for
  \code{pcgPreconditioner="nystrom"}; the default is 32 and the value is capped
  at the coefficient dimension.}
  
  \item{solveOnly}{If \code{TRUE}, no variance components are estimated: the
  mixed-model equations are solved once at known covariance parameters (genetic
  evaluation). The coefficient matrix is never assembled; a matrix-free
  preconditioned conjugate gradient iterates on the data with block-Jacobi
  preconditioning (one block for the fixed effects when there are at most 1000
  of them, and one covariance-dimension block per random-effect level).
  \code{pcgTol} is the relative-residual tolerance and \code{pcgMaxIters} the
  iteration limit (0 = automatic). Requires a Gaussian identity model,
  \code{computeCi=0}, no rotation, no \code{vcc}, a diagonal \code{W} if any,
  and a block-diagonal residual covariance whose blocks have at most 10000
  records. The result is an object of class \code{mmesSolve} with \code{b},
  \code{u}, \code{bu}, \code{uList}, \code{fitted}, \code{residuals},
  \code{theta} (covariance matrices used), \code{vcParams} (the natural-scale
  parameters used, reusable as \code{covPar}) and \code{pcg} diagnostics. No
  PEV or likelihood is computed. The default is \code{FALSE}. The fitting
  controls \code{nIters}, \code{REML}, \code{emWeight}, \code{stepWeight},
  \code{solver}, \code{henderson} and the stochastic PCG trace/log-determinant
  controls do not select or tune this separate solver.}
  
  \item{covPar}{Known covariance parameters for \code{solveOnly=TRUE}. One of:
  a fitted \code{mmes} object with the same random and residual terms (its
  final estimates are used exactly); a list with one element per term, named by
  term label or in formula order (random terms, then residual), each holding
  either the natural-scale parameter vector (as in \code{fit$covPar}) or, for
  \code{ism()}, \code{dsm()} and \code{usm()} terms, the full covariance
  matrix (e.g. \eqn{G_0} and \eqn{R_0}); or a data frame with columns
  \code{term}, \code{parameter} and \code{value} (as \code{vcParams} of an
  \code{mmesSolve} result, or the \code{vcParams} table of
  \code{returnParam=TRUE} with a \code{value} column). When \code{NULL}, every
  parameter must already be fixed in the formula (e.g.
  \code{vsm(ism(id), Gu=Ainv, sigma2=0.3, fixedSigma2=TRUE)}).}
  
}

\details{

\strong{Choosing a fitting mode}

The default, \code{henderson=TRUE}, estimates covariance parameters with
mixed-model equations in coefficient space. Supply relationship precisions
as \code{Gu}, with \code{attr(Gu, "inverse")=TRUE}. Set
\code{henderson=FALSE} to estimate parameters with the direct observation-space
engine; \code{Gu} is still supplied as a relationship precision. The direct
engine does not use the Henderson \code{solver} or PCG controls.
Non-Gaussian families use an outer PQL iteration around Gaussian working fits,
not exact GLMM maximum likelihood.

For known-parameter prediction, use \code{solveOnly=TRUE} and \code{covPar}.
This mode always uses its own matrix-free PCG solver, irrespective of
\code{solver} and \code{henderson}; it does not estimate covariance parameters
or return likelihoods or PEVs. \code{returnParam=TRUE}, by contrast, prepares
the design and covariance descriptors without fitting or solving the model.

\strong{Starting values and convergence}

When \code{vsm()} variance scales are omitted, initial scales are chosen from
fixed-effect residual variation when possible. Explicit \code{sigma2} values
are retained; \code{fixedSigma2=TRUE} prevents estimation of that scale.
Covariance constructors control their own starting shape parameters and which
of those parameters are fixed.

Henderson fits stop when either the absolute change in log-likelihood is below
\code{tolParConvLL} or the parameter convergence norm is below
\code{tolParConvNorm}. Reaching \code{nIters} alone is not convergence; inspect
\code{convergence}, \code{llik} and \code{normMonitor}. A nearly unchanged
likelihood need not imply equally precise covariance estimates. When comparing
equivalent parameterizations, tighten both convergence tolerances and allow
enough iterations. PCG trace and log-determinant approximations introduce an
additional accuracy tradeoff controlled by the PCG arguments.

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
Please review the covariance-structure vignette and pay attention to details
like format of your random and fixed variables (e.g. character and factor 
variables have different properties when returning BLUEs or BLUPs).

\strong{For tutorials} on how to perform different analysis with sommer please
look at the vignettes by typing in the terminal:

\code{vignette("sommer.covariance.structures", package="sommer")}

\code{vignette("sommer.genetic.evaluation", package="sommer")}

\code{vignette("sommer.qg", package="sommer")}

\code{vignette("sommer.gxe", package="sommer")}

\strong{Citation}

Type \emph{citation("sommer")} to know how to cite the sommer package in your publications.

\strong{Special variance structures}

 \code{\link{vsm}(\link{atm}(x,levels),ism(y))}

 can be used to specify heterogeneous variance for the "y" covariate at specific
 levels of the covariate "x", e.g., \emph{random=~vsm(atm(Location,c("A","B")),ism(ID))}
 fits a variance component for ID at levels A and B of the covariate Location.

 \code{\link{vsm}(\link{dsm}(x),ism(y))}

 can be used to specify a diagonal covariance structure for the "y" covariate 
 for all levels of the covariate "x", e.g., \emph{random=~vsm(dsm(Location),ism(ID))} 
 fits a variance component for ID at all levels of the covariate Location.

 \code{\link{vsm}(\link{usm}(x),ism(y))}

 can be used to specify an unstructured covariance structure for the "y" covariate 
 for all levels of the  covariate "x", e.g., \emph{random=~vsm(usm(Location),ism(ID))}  
 fits variance and covariance components for ID at all levels of the covariate Location.
 
 \code{\link{vsm}(\link{rrm}(Location, 2), ism(ID))}

 fits the reduced-rank covariance model with two factors across locations for
 \code{ID}. Use \code{\link{fam}(Location, 2)} instead to estimate
 location-specific variances as well as factor loadings.

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

  \code{\link{spl2Dc}(x.coord, y.coord, at.var, at.levels)}

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
https://stackoverflow.com or https://stats.stackexchange.com for the community
to help you.

If you have any bug report please go to https://github.com/covaruber/sommer or 
send me an email to address it asap, just make sure you have read the vignettes 
carefully before sending your question. If you need professional consulting 
support visit https://covaruber.github.io/sommer/

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
A Gaussian or PQL fit returns an S3 object of class \code{mmes}. The main
components are listed below; some are engine-specific or depend on
\code{computeCi}. PQL fits additionally return \code{family},
\code{dispersion}, \code{deviance} and \code{pqlMonitor}; their likelihood
and information criteria are not GLMM likelihood quantities.
With \code{returnParam=TRUE}, a setup list is returned instead, including
design matrices, covariance descriptors, \code{vcParams} and optimizer controls.
With \code{solveOnly=TRUE}, the return value has class \code{mmesSolve}; see
the \code{solveOnly} argument for its components and limitations.

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
\item{C}{the mixed-model coefficient matrix, in fixed-then-random coefficient
order. Explicit factor-score augmentation uses the augmented representation.}
\item{Ci}{the full inverse of the coefficient matrix only when
\code{computeCi=2}; otherwise it is not materialized.}
\item{theta}{a list of estimated variance covariance matrices. Each element of 
             the list corresponds to the different random and residual components}
\item{covPar}{the fitted covariance parameters in their descriptor-defined
reported coordinates.}
\item{pcgMatrixFree}{logical indicating whether the matrix-free Kronecker PCG
path was used during optimization.}
\item{residualCovarianceCompact}{logical indicating that a complete, unweighted,
non-diagonal Kronecker residual covariance was evaluated directly from its
factors without materializing the full covariance matrix. In this case the
residual element of \code{theta} is a \eqn{0 \\times 0} matrix; use
\code{covparams_mmes()} for fitted residual covariance parameters.}
\item{covParWorking}{covariance parameters on the internal optimizer scale,
such as log scales and transformed correlations. Prefer the native parameter
table for interpretation.}
\item{vcParams}{parameter identifiers and metadata for covariance descriptors,
used to specify \code{vcc} constraints.}
\item{covParNative}{a data frame containing all fixed and estimated covariance
parameters, delta-method standard errors and Z ratios in the native
model-facing scale defined by each CovarianceFactor.
For example, \code{dsm()} entries are level-specific variances rather than
variance ratios. Use \code{\link{covparams_mmes}} for estimates only or
\code{\link{covparams_mmes_se}} for estimates with standard errors.}
\item{theta_se}{estimated covariance matrix of the reported \code{covPar}
coordinates, not a vector of standard errors. Native-scale standard errors
are obtained by the delta method with \code{\link{covparams_mmes_se}}.}
\item{InfMat}{information matrix.}
\item{monitor}{covariance parameters in the reported \code{covPar}
coordinates across iterations. Explicit factor-score augmentation may instead
retain restored working-scale rows; use \code{covParNative} for interpretable
final estimates in all modes.}
\item{engine, solver}{the fitting engine and resolved Henderson solver.}
\item{engineDiagnostics}{engine-specific counters and active numerical paths;
see the \code{solver} and \code{acceleration} arguments.}
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
\item{Cchol}{legacy placeholder; the current Henderson engine returns an
empty matrix rather than exposing its internal factorization.}
\item{y}{the response vector.}
\item{W}{the column binded matrix W = [X Z] }
\item{uList}{a list containing the BLUPs in data frame format where rows are 
             levels of the random effects and column the different factors at 
             which the random effect is fitted. This is specially useful for 
             diagonal and unstructured models.}
\item{uPevList}{prediction error variances for the corresponding entries of
\code{uList}, organized by random term and covariance level when
\code{getPEV=TRUE} and \code{computeCi>0}. These are variances, not BLUPs or
standard errors.}
\item{args}{the fixed, random and residual formulas from the mmes model.}

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

## Known variances: solve the MME only (e.g. genetic evaluation)
sol1 <- mmes(Yield~Env,
             random= ~ Name + Env:Name,
             rcov= ~ units,
             data=DT, solveOnly=TRUE, covPar=ans1)
head(sol1$uList[[1]])

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
