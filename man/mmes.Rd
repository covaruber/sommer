\name{mmes}
\alias{mmes}
\title{

\strong{m}ixed \strong{m}odel \strong{e}quations \strong{s}olver

}
\description{

Fits a generalized linear mixed model using model formulas, covariance
structures and observation data. Special covariance structures can be
specified with \code{\link{vsm}()}. Numerical fitting choices are supplied
through the strictly validated \code{control} object; see
\code{\link{mmes_control}} and the \code{sommer.mmes.options} vignette.


}

\usage{

mmes(fixed, random, rcov, data, W, weights=NULL,
     naMethodX="exclude", naMethodY="exclude",
     naMethodRandom="exclude", naMethodR="exclude",
     contrasts=NULL, family=stats::gaussian(), vcc=NULL,
     covPar=NULL, control=mmes_control(), .internal=NULL, ...)
}
\arguments{
  \item{fixed}{A formula specifying the \strong{response variable(s)} \strong{and fixed effects}, i.e:

  \emph{response ~ covariate}

  By default, the fixed-effect design must have full column rank after observation
  filtering. Before fitting, \code{mmes()} reports all-zero columns (for
  example, unused levels or unobserved interaction combinations) or linearly
  dependent columns (for example, a parent factor together with unrestricted
  nested cell effects). Aliased coefficients are not automatically dropped
  or constrained to zero unless explicitly requested with
  \code{control=mmes_control(representation=list(fixedRank="drop"))}.
  This option selects a full-rank subset using shared native sparse QR for
  the Henderson and direct engines. The original coefficient layout is
  preserved, with omitted coefficients returned as zero identification
  constraints and flagged in \code{summary(fit)$betas}; these zeros are not
  estimates of independently identifiable effects. Details, mappings and
  supported inference workflows are documented in \code{\link{mmes_control}}.

  Apply \code{droplevels()} to the retained data to remove unused factor
  levels. For a sparse interaction, an observed-cell factor created with
  \code{interaction(..., drop=TRUE)} and a formula such as
  \code{response ~ 0 + cell} can represent the same fitted cell means without
  redundant columns. Adding the parent factor again reintroduces aliasing.
  Dropping unused levels alone does not remove all unobserved interaction
  combinations. The diagnostic setup mode
  (\code{control=mmes_control(output=list(setup=TRUE))}) still returns the
  design without requiring full column rank, allowing it to be inspected.
  With \code{fixedRank="drop"}, setup instead returns the reduced design and
  its \code{fixedRank} mapping.

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

  \code{\link{rrm}} for a reduced-rank random covariance with specific
  variances fixed at zero; use \code{fam(..., variance="homogeneous")}
  for an estimated common specific variance;

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
  \code{vsm(..., Gu=Gu, rotation=TRUE)}. With
  \code{control$formulation="henderson"}, random
  coefficients are represented in the eigenbasis and the diagonal relationship
  precision is used in the mixed model equations. With
  \code{control$formulation="direct"},
  the equivalent transformed observation-covariance parameterization is used.
  Public random effects, fitted values, and residuals are returned in the
  original basis. The current implementation requires a complete balanced
  Gaussian layout, one rotated relationship term, and
  \code{control$output$inverse=0} during
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
  \item{W}{Weights matrix (e.g., when covariance among plots exist). Internally
  W is squared and inverted as Wsi = solve(chol(W)), then the residual matrix
  is calculated as R = Wsi*O*Wsi.t(), where * is the matrix product, and O 
  is the original residual matrix. Weights can be combined with any residual
  structure, including Kronecker (multi-trait) residuals.}
  \item{weights}{Optional one-sided formula describing independent row blocks
  of the supplied \code{W} matrix, for example \code{weights = ~trial}. Rows
  with different formula groups must have zero cross-group entries in \code{W}.
  This is structural metadata only: it does not generate numeric weights.
  Declared blocks are processed independently by the Henderson weighting
  path; without this formula, \code{W} is handled without a known grouping.}
     \item{control}{A \code{sommer_mmes_control} object created by
  \code{\link{mmes_control}}. It owns numerical task, equation formulation,
  representation, solver, determinant, score, curvature, covariance update,
  proposal, convergence, output, execution and diagnostic settings. Its default
  is \code{mmes_control()}. Removed numerical arguments are not accepted by
  \code{mmes()}; configure those choices in this object. See the control help
  and vignette for defaults, mathematical ownership, eligibility and migration.}

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

   \item{contrasts}{an optional list. See the contrasts.arg of model.matrix.default.}
  \item{vcc}{Optional data frame of equality/scaling constraints between
  covariance parameters (columns \code{parameter}, \code{group} and optionally
  \code{scale}), identified through the \code{vcParams} table returned by
  \code{control$output$setup=TRUE} and stored in the fit. See \code{\link{vcc}}.}

  \item{family}{A \code{stats::family()} object. The default Gaussian identity
  family uses the ordinary linear mixed-model implementation. Other families are
  fitted by penalized quasi-likelihood (PQL): each outer iteration forms an IRLS
  working response and weights, then fits the resulting weighted Gaussian mixed
  model with the Henderson solver. The initial implementation supports one numeric
  response. Any family with a valid link can be used, e.g. \code{binomial()}
  (logit, probit, cloglog), \code{poisson()}, \code{Gamma()},
  \code{inverse.gaussian()}, \code{gaussian(link="log")}, \code{quasipoisson()}
  and \code{MASS::negative.binomial(theta)}. Binomial counts out of a known
  number of trials are supplied as a single proportion response built with
  \code{\link{binm}}, whose trials are used as prior weights. Binomial,
  Poisson and negative binomial working residual dispersions are fixed to one.
  PQL uses outer-loop controls from \code{control$task$pql}; compatibility with
  update and solver requests is model-dependent.
  For other families the reported \code{dispersion} is the Pearson estimate.
  The negative binomial \eqn{\theta} is estimated (conditional on the random
  effects) with \code{control$task$pql=list(estimateTheta=TRUE)} and returned in
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

    \item{covPar}{Known covariance parameters for the solve task selected by
  \code{control$task$mode="solve"}. Accepted forms are a compatible fitted
  \code{mmes} object, a named list of term parameters, or a parameter data
  frame. When \code{NULL}, covariance values must be fixed in the model.}

  \item{.internal}{Reserved for package-internal recursive fitting state.
  Users should not set this argument.}

 \item{...}{Deprecated legacy arguments. The supported arguments are
 \code{getPEV}, \code{tolParInv}, and \code{verbose}; they map respectively to
 \code{control$output$inverse}, \code{control$execution$inverseTol}, and
 \code{control$diagnostics$verbose}. These arguments override matching
 settings in \code{control} and issue a deprecation warning.}
}
\details{

\strong{Structurally unestimable REML components}

For Gaussian identity-link REML fits, \code{mmes()} checks whether an
independently scaled random component contributes anything outside the
fixed-effect design space, using the observations retained by the common
missing-data filter. A component whose entire random design lies in that
space is invisible to the restricted likelihood. Such components are
removed before optimization, with a warning, rather than assigned an
arbitrary variance by a regularized information-matrix solve.

The conservative check handles a whole random term with a free variance
scale and fixed covariance shape, or individual levels of a single
\code{dsm()} or \code{atm()} shaping factor when its scale and all variance
ratios are free. For example, with fixed locale means, a locale-specific
independent block variance is omitted when that locale has only one
observed block. A shared block variance is retained if other locales
provide information about it. Correlated, multiply shaped, constrained or
partly fixed structures are not automatically reduced at individual levels.
The check does not diagnose general aliasing between covariance terms or
weak identification near a variance boundary.

Automatic omission is not applied to ML, PQL, known-covariance solve tasks,
\code{vcc} constraints, observation rotation or factor-score augmentation.
For supported fits, omitted components retain their original positions in
random-effect and covariance outputs, with zero covariance and zero random
effects. These zeros denote an omission convention, not variances estimated
to be zero. Their native-scale standard errors and Z-ratios are \code{NA};
\code{summary()} and \code{covparams_mmes()} label them
\code{"unestimable (omitted)"}. Only retained free parameters contribute to
optimization and information criteria. Omission preserves the restricted
likelihood information about retained parameters; it does not establish
zero between-block variation in the population. Inference on fixed means
and predictions is conditional on the omission convention.

\strong{Fixed factors with many levels}

In Henderson fits the fixed-effect design \code{X} is stored as a sparse
matrix. The coefficient matrix \code{C} of the mixed model equations is
factorized with sparse methods at every REML iteration, whatever the
\code{control$output$inverse} setting. The fill-reducing ordering of the
sparse factorization eliminates the fixed-effect equations within the same
factor as the random ones, which is equivalent to absorbing them. A fixed
factor with thousands of levels, such as contemporary groups, does not
create a dense fixed-effect block during estimation. Variance components,
BLUEs and BLUPs are identical to those of an explicit absorption of the
fixed effects.

\code{control$output$inverse} only controls work done after convergence:
\describe{
  \item{\code{0} (default)}{No inverse of \code{C} is formed. Fixed-effect
  standard errors in \code{summary()} are \code{NA}.}
  \item{\code{1}}{Elements of the inverse are computed only within the
  sparsity pattern of the factor (Takahashi selected inverse), giving
  random-effect PEVs. The fixed part stays sparse. Fixed-effect standard
  errors in \code{summary()} remain \code{NA}.}
  \item{\code{2}}{The complete inverse of \code{C}, including the fixed
  block, is formed as a dense matrix. This needs about
  \eqn{8(p+q)^2} bytes for \eqn{p} fixed and \eqn{q} random effects
  (around 1.4 GB for 13,000 effects) and is the mode that reports
  fixed-effect standard errors in \code{summary()}.}
}
With many fixed levels, keep \code{inverse=0} or \code{1}. When standard
errors of particular fixed-effect combinations are needed, use
\code{\link{predict.mmes}}, which computes exact standard errors from the
sparse \code{C} without forming \code{Ci}. \code{\link{postPEV}} can add
PEVs or the full inverse after fitting.

When even the sparse REML fit is too slow, a practical alternative is to
estimate the variance components with the large factor fitted as random,
optionally on a well-connected subset of the data. Then refit with the
factor as fixed using \code{control=mmes_control(task=list(mode="solve"))}
and \code{covPar} to obtain BLUEs and BLUPs with a single factorization.
Simply omitting the factor during REML is not equivalent: its variation is
then attributed to the residual and to confounded random terms.

\strong{Numerical controls}

The opt-in \code{control$representation$equations="conditional"} path fits Gaussian
FA models using sparse latent equations with exact REML AI updates by default,
or conditional-moment EM through \code{control$update$method="em"}. The earlier
\code{control$update$method="conditional-fa"} selector retains its EM behavior.
Both optimizers support incomplete and
unbalanced observations. It accepts one FA relationship term and diagonal
or FA residuals. Matching \code{fam(shared="group")} tags request common
trait loadings across genetic and residual sources. The default AI path is
unchanged. See \code{\link{mmes_control}} and the conditional FA section of
\code{vignette("sommer.mmes.options")} for eligibility and output restrictions.
Conditional fits retain factorized covariance descriptors, not dense
\code{theta} matrices, public \code{C}, coefficient uncertainty, or
covariance-parameter uncertainty. Use \code{covmatrix_mmes(..., se=FALSE)}
to explicitly reconstruct a bounded-size covariance.

\code{mmes()} keeps model specification separate from numerical policy.
The \code{control} argument defaults to \code{mmes_control()}, whose
modules have strict schemas: unknown field names are errors. The default
configuration retains the former standard bundle: Henderson equations, REML,
exact AI curvature, EM recovery, no extrapolation, 30 inner iterations,
full covariance information, coefficient-matrix retention, and inherited
solver/resource defaults. The \code{mmes_control} help gives every field,
default, module owner and validation rule.

\strong{Ownership and eligibility}

\describe{
  \item{Task}{\code{control$task} selects fitting or known-covariance solving,
  REML or ML, and outer PQL controls.}
  \item{Formulation}{\code{control$formulation} selects Henderson or direct
  observation-space equations.}
  \item{Representation and solve}{\code{control$representation} describes
  equation/covariance representation; \code{control$solve} owns linear solver,
  preconditioning, reuse and elimination requests.}
  \item{Evaluation}{\code{control$logdet} owns determinant methods;
  \code{control$score} owns score traces and covariance checks;
  \code{control$curvature} owns exact or hybrid AI and reduced-direction work.}
  \item{Optimization}{\code{control$update} owns EM/PXEM recovery, blending
  and common-scale profiling. \code{control$proposal} owns step scaling and
  Aitken extrapolation. Aitken and hybrid AI are mutually exclusive.}
  \item{Convergence and outputs}{\code{control$convergence} owns inner
  iteration, tolerances and stochastic accuracy; \code{control$output} owns
  coefficient/covariance uncertainty and setup output.}
  \item{Execution and diagnostics}{\code{control$execution} owns memory,
  hardware, threading, precision and seeds. \code{control$diagnostics} owns
  verbosity, package-age warnings and profiling.}
}

These modules do not form a Cartesian product of interchangeable algorithms.
Model-dependent admission remains authoritative; a requested mode may be
ineligible. Consult \code{\link{mmes_control}} for prerequisites, exclusions
and current behavior. No new method combinations are implied by the modular API.
The internal \code{.mmes_capabilities()} catalog is descriptive and does
not resolve eligibility or predict execution.

}

\value{
If all parameters are correctly indicated the program will return a list with the following information:

\item{data}{the dataset used in model fitting after application of the common
              observation mask.}
\item{obsInfo}{an observation-level audit table describing the common missing-data
               filter. It includes the original row index, validity indicators for
               response, fixed, random, and residual model components, the final
               inclusion indicator, and the reason for exclusion when applicable.}
\item{remlOmissions}{when present, a data frame identifying structurally
unestimable random covariance components omitted from REML, with columns
\code{term}, \code{parameter} and \code{reason}. It is also available in
\code{engineDiagnostics$remlOmissions}.}
\item{Dtable}{the table to be used for the predict function to help the program
              recognize the factors available.}
\item{llik}{the vector of log-likelihoods across iterations}
\item{b}{the vector of fixed effects on the original response scale.}
\item{u}{the vector of random effects on the original response scale, consistent
with \code{uList}.}
\item{bu}{the vector of fixed and random effects together on the original response
scale, formed by concatenating \code{b} and \code{u}.}
\item{Ci}{the inverse of the coefficient matrix.}
\item{Ci_11}{the inverse of the coefficient matrix pertaining to the fixed effects.}
\item{theta}{a list of estimated variance covariance matrices. Each element of 
             the list corresponds to the different random and residual components}
\item{covPar}{the fitted covariance parameters in their descriptor-defined
reported coordinates.}
\item{control}{the validated modular numerical controls used for the fit.
See \code{\link{mmes_control}}.}
\item{pcgRecycleState}{present only when \code{pcg$exportRecycle=TRUE}.
A versioned serializable Krylov basis with coefficient-coordinate metadata,
for explicit import into a compatible PCG fit.}
\item{covarianceInformation}{\code{"none"} when covariance information was
explicitly omitted. The default retains full covariance information.}
\item{retainCoefficientMatrix}{FALSE when an opt-in C-free saddle-point fit
omitted the final compatibility matrix. In that mode \code{C} is NULL, and
C-based post-fit methods require a refit retaining C. The default retains C.}
\item{engineDiagnostics}{path-dependent numerical diagnostics.
\code{compression} reports method, original observations, groups, residual
classes, retained statistic bytes and setup time. \code{pcg} reports iterative
work, reuse and stochastic accuracy diagnostics for REML PCG.
\code{profile} contains detailed phase timings when explicitly enabled.
Retained statistic bytes are not process peak memory.}
\item{pcgMatrixFree}{logical indicating whether the matrix-free Kronecker PCG
path was used during optimization.}
\item{residualCovarianceCompact}{logical indicating that a complete, unweighted,
non-diagonal Kronecker residual covariance was evaluated directly from its
factors without materializing the full covariance matrix. In this case the
residual element of \code{theta} is a \eqn{0 \\times 0} matrix; use
\code{covparams_mmes()} for fitted residual covariance parameters.}
\item{covParNative}{a data frame containing all fixed and estimated covariance
parameters in the native model-facing scale defined by each CovarianceFactor.
For example, \code{dsm()} entries are level-specific variances rather than
variance ratios. The table includes delta-method \code{StdError} and
\code{Zratio} columns when covariance information is retained, matching
\code{covparams_mmes_se()}. With \code{covarianceInformation="none"} these
columns are NA; \code{covparams_mmes()} still returns point estimates.}
\item{theta_se}{covariance-parameter uncertainty matrix; NULL when
\code{covarianceInformation="none"}.}
\item{InfMat}{information matrix; NULL when
\code{covarianceInformation="none"}.}
\item{monitor}{The values of the variance-covariance components across iterations 
              during the REML estimation.}
\item{AIC}{Akaike information criterion}
\item{BIC}{Bayesian information criterion}
\item{convergence}{a TRUE/FALSE statement indicating if the model converged.}
\item{partitions}{a list where each element contains a matrix indicating where 
                  each random effect starts and ends.}
\item{partitionsX}{a list where each element contains a matrix indicating where 
                   each fixed effect starts and ends.}
\item{percDelta}{the matrix of percentage change in deltas (see \code{control$convergence$parameterTol}).}
\item{normMonitor}{the matrix of the three norms calculated (see \code{control$convergence$parameterTol}).}
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

## Known variances: solve the MME only (e.g. genetic evaluation)
sol1 <- mmes(Yield~Env,
             random= ~ Name + Env:Name,
             rcov= ~ units,
             data=DT, control=mmes_control(task=list(mode="solve")),
             covPar=ans1)
head(sol1$uList[[1]])

\donttest{

####===========================================####
#### Univariate heterogeneous variance models  ####
####===========================================####
DT=DT[with(DT, order(Env)), ]
## Compound simmetry (CS) + Diagonal (DIAG) model
ans2 <- mmes(Yield~Env,
             random= ~Name + vsm(dsm(Env),ism(Name)),
             rcov= ~ dsumm(vsm(ism(units)), by=Env),
             data=DT)
summary(ans2)

####===========================================####
####  Univariate unstructured variance models  ####
####===========================================####

ans3 <- mmes(Yield~Env,
             random=~ vsm(usm(Env),ism(Name)),
             rcov= ~ dsumm(vsm(ism(units)), by=Env),
             data=DT)
summary(ans3)


}

\donttest{

## simulated MET
nInds=300
nMarks=1000
nEnvs=3

#random population of nInds lines with nMarks markers
M <- matrix(rep(0,nInds*nMarks),nInds,nMarks)
for (i in 1:nInds) {
  M[i,] <- ifelse(runif(nMarks)<0.5,-1,1)
}
# compute the relationship matrix
MMT <- tcrossprod(M)
m <- sum(diag(MMT))/nrow(MMT)
MMT <- MMT/m
MMT <- MMT + diag(1e-05, ncol(MMT), ncol(MMT))
# get phenotypes
rg <- enhancer::simGECorMat(nEnvs, 1, mu=0.2, v=0.2)
rg <- Matrix::nearPD(rg)$mat
rg # desired genetic correlation
u <- MASS::mvrnorm(n = nMarks, rep(0, nEnvs), Sigma=rg) # marker effects
g <- crossprod(t(M),u) # genetic values
h2 <- 0.2  #heritability 
envEffect <- rnorm(nEnvs, mean = abs(rnorm(1))) # environment means
e <- MASS::mvrnorm(n = nInds, envEffect, Sigma=(1-h2)/h2*diag(diag(var(g)))) # errors
p = g + e # phenotypes
y <- as.vector(p) # phenotypes as vector
env <- sort(rep(paste0("e",1:nEnvs), nrow(g))) # environment labels
ind <- rep(paste0("i",1:nrow(g)), nEnvs ) # individuals labels
colnames(MMT) <- rownames(MMT) <- paste0("i",1:nrow(g)) # assign dimnames

MMTi <- solve(MMT + diag(1e-4,ncol(MMT),ncol(MMT)))
MMTi[lower.tri(MMTi)] <- t(MMTi)[lower.tri(MMTi)] # fill the lower triangular
MMTi <- as(as(as( MMTi,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
attr(MMTi, 'inverse')=TRUE

DT <- data.frame(ind=ind, env=env, y=y)
# fit the models
mixm2 <- mmes(y ~ env,
              random=~ vsm(dsm(env), ism(ind), Gu=MMTi),
              data=DT
)
summary(mixm2)$varcomp

}

}
