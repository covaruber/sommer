# **m**ixed **m**odel **e**quations **s**olver

Fits a generalized linear mixed model coded in C++ using the Armadillo
and Eigen libraries to optimize matrix operations. Special covariance
structures can be specified with the
[`vsm()`](https://covaruber.github.io/sommer/reference/vsm.md) function.
Two engines are available; henderson formulation and direct reml
(henderson argument) and three solvers in each; ldlt, cholmod and pcg
(solver argument) to adjust to the modeling situation.

## Usage

``` r
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
    factorScoreAugmentation="none", .factorScoreParameters=NULL,
    pcgPreconditioner="diagonal", pcgNystromRank=32L,
    solveOnly=FALSE, covPar=NULL)
```

## Arguments

- fixed:

  A formula specifying the **response variable(s)** **and fixed
  effects**, i.e:

  *response ~ covariate*

- random:

  A formula specifying the **random effects**, e.g., *random = ~
  genotype + year*. A simple random term that is not wrapped in
  [`vsm`](https://covaruber.github.io/sommer/reference/vsm.md) is
  treated as an identity-covariance random effect.

  The [`vsm`](https://covaruber.github.io/sommer/reference/vsm.md)
  function is the main interface for specifying structured covariance
  models for random effects. Its current covariance model is

  \$\$G = \sigma^2 (K_1 \otimes K_2 \otimes \cdots \otimes K_m) \otimes
  A,\$\$

  where \\\sigma^2\\ is the single variance scale owned by
  [`vsm()`](https://covaruber.github.io/sommer/reference/vsm.md),
  \\K_1,\ldots,K_m\\ are covariance-shaping factors supplied by
  covariance constructors, and \\A\\ is the known covariance
  relationship for the final random effect. In the Henderson
  implementation `Gu` supplies the **inverse** of \\A\\; it must have
  row and column names matching the levels of the final random effect
  and `attr(Gu,"inverse")=TRUE`.

  Within `vsm(...)`, the **last** constructor supplies the incidence
  matrix for the main random effect, while all preceding constructors
  define covariance-shaping factors. The number of covariance factors is
  not limited. For example,

  `random = ~ vsm(dsm(Location), ism(Name))`

  fits heterogeneous variances across `Location` for the random effect
  `Name`, while

  `random = ~ vsm(dsm(Environment), ar1m(Row), ar1m(Column), ism(Name), Gu=Ainv)`

  specifies an arbitrary Kronecker product of covariance structures
  before the relationship matrix for `Name`.

  Covariance constructors currently available include:

  [`ism`](https://covaruber.github.io/sommer/reference/ism.md) for
  identity covariance;

  [`dsm`](https://covaruber.github.io/sommer/reference/dsm.md) and
  [`atm`](https://covaruber.github.io/sommer/reference/atm.md) for
  diagonal and selected-level diagonal covariance structures;

  [`usm`](https://covaruber.github.io/sommer/reference/usm.md) for an
  unstructured covariance matrix;

  [`csm`](https://covaruber.github.io/sommer/reference/csm.md) and
  [`corgm`](https://covaruber.github.io/sommer/reference/corgm.md) for
  correlation structures;

  [`ar1m`](https://covaruber.github.io/sommer/reference/ar1m.md),
  [`ar2m`](https://covaruber.github.io/sommer/reference/ar2m.md), and
  [`ar3m`](https://covaruber.github.io/sommer/reference/ar3m.md) for
  autoregressive covariance structures;

  [`mam`](https://covaruber.github.io/sommer/reference/mam.md)
  (including `ma1m` and `ma2m`) for moving-average covariance
  structures;

  [`toeplitzm`](https://covaruber.github.io/sommer/reference/toeplitzm.md)
  for a general positive-definite Toeplitz correlation structure;

  [`fam`](https://covaruber.github.io/sommer/reference/fam.md) for
  factor-analytic covariance;

  [`rrcm`](https://covaruber.github.io/sommer/reference/rrcm.md) for the
  reduced-rank covariance approximation;

  [`antem`](https://covaruber.github.io/sommer/reference/antem.md) for
  antedependence covariance;

  [`maternm`](https://covaruber.github.io/sommer/reference/maternm.md)
  for Matern spatial covariance;

  [`sar`](https://covaruber.github.io/sommer/reference/sar.md) and
  [`car`](https://covaruber.github.io/sommer/reference/car.md) for
  simultaneous and conditional autoregressive spatial covariance models;
  and

  [`ownm`](https://covaruber.github.io/sommer/reference/ownm.md) for a
  fixed or user-defined covariance function.

  Design-matrix utilities such as
  [`overlay`](https://rdrr.io/pkg/enhancer/man/overlay.html),
  [`spl2Dc`](https://covaruber.github.io/sommer/reference/spl2Dc.md),
  [`leg`](https://rdrr.io/pkg/enhancer/man/leg.html), and
  [`redmm`](https://rdrr.io/pkg/enhancer/man/redmm.html) can still be
  combined with
  [`vsm()`](https://covaruber.github.io/sommer/reference/vsm.md) when
  appropriate.

  See [`vsm`](https://covaruber.github.io/sommer/reference/vsm.md) and
  the individual covariance-constructor help pages for details on
  parameterizations and examples.

  A single relationship term may request Lee–van der Werf rotation with
  `vsm(..., Gu=Gu, rotation=TRUE)`. With `henderson=TRUE`, random
  coefficients are represented in the eigenbasis and the diagonal
  relationship precision is used in the mixed model equations. With
  `henderson=FALSE`, the equivalent transformed observation-covariance
  parameterization is used. Public random effects, fitted values, and
  residuals are returned in the original basis. The current
  implementation requires a complete balanced Gaussian layout, identity
  residual covariance, default identity `W`, one rotated relationship
  term, and `computeCi=0` during fitting.

- rcov:

  A formula specifying the **residual covariance structure**. The
  default is *rcov = ~ units*, which fits a homogeneous residual
  variance.

  Structured residual covariance is specified with the same
  [`vsm`](https://covaruber.github.io/sommer/reference/vsm.md) and
  CovarianceFactor constructors used for random effects. The residual
  term should be expressed as a single
  [`vsm()`](https://covaruber.github.io/sommer/reference/vsm.md) term;
  arbitrary covariance complexity is represented by placing multiple
  covariance factors inside that term rather than by summing multiple
  residual terms.

  The final term is normally `ism(units)`, where `units` is the
  observation-level incidence generated internally by `mmes`. For
  example,

  `rcov = ~ vsm(dsm(Environment), ism(units))`

  fits heterogeneous residual variances across environments, while

  `rcov = ~ vsm(ar1m(Row), ar1m(Column), ism(units))`

  fits a separable row-by-column AR1 residual covariance. More
  generally, structures such as
  [`usm`](https://covaruber.github.io/sommer/reference/usm.md),
  [`csm`](https://covaruber.github.io/sommer/reference/csm.md),
  [`toeplitzm`](https://covaruber.github.io/sommer/reference/toeplitzm.md),
  [`maternm`](https://covaruber.github.io/sommer/reference/maternm.md),
  [`sar`](https://covaruber.github.io/sommer/reference/sar.md),
  [`car`](https://covaruber.github.io/sommer/reference/car.md), and
  [`ownm`](https://covaruber.github.io/sommer/reference/ownm.md) can be
  used when their incidence layout is appropriate for the residual
  model.

  Residual covariance factors must identify exactly one
  covariance-product coordinate for each observation. `mmes()`
  constructs the residual block and local-coordinate indexing
  internally, so the data do **not** need to be manually sorted by the
  variables defining the residual structure.

  As with random effects,
  [`vsm()`](https://covaruber.github.io/sommer/reference/vsm.md) owns
  one residual variance scale \\\sigma^2\\; the preceding covariance
  constructors provide dimensionless covariance shapes that are combined
  through an arbitrary Kronecker product.

  See [`vsm`](https://covaruber.github.io/sommer/reference/vsm.md) and
  the individual covariance-constructor help pages for details.

- data:

  Optional data frame containing variables used by the model. Model
  expressions are evaluated using normal R scoping rules: variables are
  first sought in `data` and may also be obtained from the environments
  associated with the model formulas or from the calling environment.
  Consequently, auxiliary objects such as relationship matrices,
  coordinate vectors, or grouping variables do not have to be copied
  into `data`. If `data` is omitted, variables are evaluated from the
  formula/calling environment.

  Observation-level variables used by the response, fixed effects,
  random effects, and residual covariance specification must
  nevertheless have compatible lengths. Internally, `mmes()` constructs
  a common observation mask and applies it consistently to all
  observation-level model components.

- W:

  Weights matrix (e.g., when covariance among plots exist). Internally W
  is squared and inverted as Wsi = solve(chol(W)), then the residual
  matrix is calculated as R = Wsi\*O\*Wsi.t(), where \* is the matrix
  product, and O is the original residual matrix.

- weights:

  Optional one-sided formula declaring independent row blocks of `W`,
  for example `weights=~trial`. Rows assigned to different groups must
  have zero cross-group entries in `W`; the formula supplies structural
  metadata and does not create numeric weights.

- factorScoreAugmentation:

  `"none"` uses the marginal covariance MME. `"fixed-shape"` and
  `"profile"` can use FA/RR or nonnegative compound-symmetry random
  factors. `usm()` remains on the marginal path because it is full-rank.
  Augmentation requires Henderson REML, `computeCi=0`, one eligible
  shaping factor per random term, no rotation, and no user `vcc`.

- nIters:

  Maximum number of iterations allowed.

- tolParConvLL:

  Convergence criteria based in the change of log-likelihood between
  iteration i and i-1.

- tolParConvNorm:

  When using the Henderson method this argument is the convergence
  criteria based in the norm proposed by Jensen, Madsen and Thompson
  (1997):

  e1 = \|\| InfMatInv.diag()/sqrt(N) \* dLu \|\|

  where InfMatInv.diag() is the diagonal of the inverse of the
  information matrix, N is the total number of variance components, and
  dLu is the vector of first derivatives.

- tolParInv:

  Tolerance parameter for matrix inverse used when singularities are
  encountered in the estimation procedure. By default the value is
  1e-06. This parameter should be fairly small because the it is used to
  bend matrices like the information matrix in the henderson algorithm
  or the coefficient matrix when it is not positive-definite.

- naMethodX:

  Missing-data policy for variables entering the fixed-effect design
  matrix. The default, `"exclude"`, marks observations with missing
  fixed-effect information for removal from the common observation set.
  `"include"` retains the historical sommer behavior for supported cases
  by imputing missing fixed-effect covariates. This option should be
  used carefully, because imputation changes the fitted design matrix.

- naMethodY:

  Missing-data policy for the response variable(s). The default,
  `"exclude"`, removes observations with missing response information
  from the common observation set. The historical `"include"` and
  `"include2"` behaviors are retained where supported: `"include"`
  imputes missing responses, whereas `"include2"` is intended for
  multivariate-response situations in which an observation may be
  retained when at least one response is observed.

- naMethodRandom:

  Missing-data policy for observation-level variables used to construct
  random-effect terms. The default is `"exclude"`. Under this policy, an
  observation for which a required random-effect classification or
  covariate cannot be evaluated is excluded from the common observation
  set. Random terms are evaluated before the final observation subset is
  applied, so the same inclusion mask is used consistently for the
  response, fixed-effect design, and all random-effect design matrices.

- naMethodR:

  Missing-data policy for observation-level variables used to define the
  residual covariance structure. The default is `"exclude"`.
  Observations lacking information required to assign their residual
  covariance-product coordinate are excluded from the common observation
  set. The residual block and local-coordinate layout are then
  constructed on the retained observations.

- returnParam:

  A TRUE/FALSE value to indicate if the program should return the
  parameters to be used for fitting the model instead of fitting the
  model.

- dateWarning:

  A TRUE/FALSE value to indicate if the program should warn you when is
  time to update the sommer package.

- verbose:

  A TRUE/FALSE value to indicate if the program should return the
  progress of the iterative algorithm.

- stepWeight:

  A vector of values (of length equal to the number of iterations)
  indicating the weight used to multiply the update (delta) for variance
  components at each iteration. If NULL the 1st iteration will be
  multiplied by 0.5, the 2nd by 0.7, and the rest by 0.9. This argument
  can help to avoid that variance components go outside the parameter
  space in the initial iterations which happens very often with the AI
  method but it can be detected by looking at the behavior of the
  likelihood. In that case you may want to give a smaller weight.

- emWeight:

  A vector of values (of length equal to the number of iterations)
  indicating with values between 0 and 1 the weight assigned to the EM
  information matrix; the values 1 - emWeight are applied to the AI
  information matrix to produce the joint information matrix used in the
  variance-component update. This creates an EM-heavy warm start early
  in the fit and gradually transitions to an AI-dominated update as the
  optimizer approaches convergence. The default schedule is a smooth
  decay from 1 at the first iteration to 0.05 by the last iteration,
  i.e. the fit begins in a strongly stabilized EM regime and then moves
  toward the classical AI update. Values outside `[0, 1]` are rejected.

- contrasts:

  an optional list. See the contrasts.arg of model.matrix.default.

- getPEV:

  A logical value indicating whether prediction error variance results
  should be organized and returned when they are requested through
  `computeCi`.

- henderson:

  A logical value indicating which REML engine to use. `TRUE` (the
  default) uses the Henderson mixed-model-equations algorithm, efficient
  when there are many more records than coefficients to estimate.
  `FALSE` uses a direct-inversion engine that inverts the phenotypic
  covariance matrix in observation space instead; it is less efficient
  than `henderson=TRUE` when there are many more records than
  coefficients, but can be more efficient when there are many more
  coefficients to estimate than records available (e.g.,
  marker/SNP-BLUP-style models). Both engines share the same
  [`vsm()`](https://covaruber.github.io/sommer/reference/vsm.md)
  covariance-structure parameterization. The direct-inversion engine
  currently supports `computeCi` values `0` and `2` only (no Takahashi
  selected-inverse mode), has no `solver`/PCG/CHOLMOD options, and
  requires a single response column (use the long-format
  `vsm(usm(trait), ...)` convention for multi-trait models, as with the
  Henderson engine).

- computeCi:

  An integer indicating whether post-fit prediction error variance
  information and/or the inverse of the mixed model coefficient matrix
  should be computed when the Henderson algorithm is used
  (`henderson=TRUE`). The available options are:

  `0`: no additional inverse-related computation is performed after
  convergence. The full inverse of the coefficient matrix (`Ci`) is not
  formed and `uPevList` is not computed. This is the fastest option and
  is the default.

  `1`: prediction error variances are computed using the Takahashi
  sparse inverse subset algorithm. This method uses the sparse LDLT
  factorization of the coefficient matrix and recursively evaluates only
  those elements of the inverse that belong to the sparsity pattern
  induced by the factorization. In particular, the diagonal elements
  required for prediction error variances are obtained without
  constructing the complete inverse matrix. Therefore, `uPevList` is
  returned while `Ci` is not fully materialized.

  With `solver="cholmod"` (see below), the supernodal factorization used
  during REML iterations has no Takahashi selected-inverse equivalent.
  In that case `computeCi=1` automatically performs one additional
  sparse LDLT factorization of the converged coefficient matrix after
  REML convergence, purely to extract the selected inverse subset. This
  is a single one-time cost, not repeated every iteration, so the
  iterative speed benefit of `solver="cholmod"` is retained.

  `2`: the complete inverse of the coefficient matrix is computed. The
  prediction error variances in `uPevList` are then extracted from the
  appropriate diagonal elements of `Ci`. This option is the most
  computationally and memory intensive, especially for large mixed model
  equation systems. With `solver="cholmod"` the full inverse is obtained
  directly from the supernodal factor (no extra LDLT refactorization is
  needed for this mode).

  For large models, `computeCi=0` is recommended during model fitting.
  If prediction error variances or the full inverse are needed
  afterwards, they can be obtained with
  [`postPEV`](https://covaruber.github.io/sommer/reference/postPEV.md)
  without refitting the model.

- solver:

  Linear-system solver used by the Henderson implementation. `"auto"`
  (default) picks a solver automatically based on the density of the
  random-effect relationship matrices (`Gu`) actually supplied: if any
  random effect uses a dense `Gu` (e.g. a genomic or marker-based
  relationship matrix, which is typically close to fully dense),
  `"cholmod"` is selected; otherwise (e.g. sparse pedigree-based
  relationship matrices, or no `Gu` at all) `"ldlt"` is selected.
  Passing an explicit value (`"ldlt"`, `"pcg"`, or `"cholmod"`) always
  overrides the automatic choice. `"ldlt"` uses the sparse direct
  simplicial LDLT factorization path. `"pcg"` selects the iterative
  preconditioned conjugate-gradient path. `"cholmod"` uses the
  supernodal (blocked, BLAS-3) sparse Cholesky factorization bundled
  with R's `Matrix` package (SuiteSparse CHOLMOD), accessed through
  `Matrix`'s public C API - no additional system library needs to be
  installed, since `Matrix` is already a dependency of `sommer`.
  Supernodal factorization can be faster than the simplicial LDLT path
  for large coefficient matrices (e.g., models built on large pedigree
  or genomic relationship matrices), because it reorganizes the sparse
  factorization into dense blocks that are handled with threaded BLAS-3
  kernels rather than column-by-column updates; the speed-up depends on
  the BLAS linked to the running R installation. `solver="cholmod"`
  currently supports all `computeCi` values (see above for how
  `computeCi=1` is handled for this solver). The covariance-model
  specification is independent of this choice: all solvers receive the
  same generic CovarianceFactor descriptors produced by
  [`vsm()`](https://covaruber.github.io/sommer/reference/vsm.md).

- pcgTol:

  Convergence tolerance used by the PCG linear solver. The default is
  `1e-8`. This argument is used when `solver="pcg"`. Smaller values
  request more accurate iterative solves but may increase computation
  time.

- pcgMaxIters:

  Maximum number of PCG iterations for an individual linear solve. A
  value of `0` lets the C++ solver choose its internal iteration limit
  from the problem dimension. This argument is used only by the PCG
  path.

- pcgTraceProbes:

  Number of stochastic probe vectors used by the PCG-based machinery for
  trace quantities required during REML calculations. The default is
  `8`. Increasing this value can reduce stochastic approximation
  variability at additional computational cost. This argument is
  relevant to `solver="pcg"`.

- pcgLanczosSteps:

  Number of Lanczos steps used by the PCG-based stochastic
  log-determinant calculations. The default is `20`. Larger values can
  improve the approximation at additional computational cost. This
  argument is relevant to `solver="pcg"`.

- REML:

  Logical, `TRUE` (default) fits variance/covariance parameters by
  restricted maximum likelihood (REML), the classical choice for
  BLUP-style variance component estimation. Setting `REML=FALSE` fits by
  maximum likelihood (ML) instead: the log-likelihood, score equations,
  and reported `monitor`/log-lik values drop the REML "restriction" term
  that accounts for the degrees of freedom used to estimate the fixed
  effects. ML estimates of variance components are known to be biased
  downward relative to REML, but unlike REML log-likelihoods, ML
  log-likelihoods are directly comparable (via a likelihood-ratio test)
  across models that differ in their **fixed**-effects specification -
  REML likelihoods are only comparable across models that share the same
  fixed effects. `REML=FALSE` currently requires `solver="ldlt"` or
  `solver="cholmod"` (`solver= "auto"` already resolves to one of
  these). `solver="pcg"` is not yet supported with `REML=FALSE`.
  Prediction error variances (`uPevList`/`computeCi`) are still computed
  from the REML-style inverse regardless of `REML`.

- family:

  A [`stats::family()`](https://rdrr.io/r/stats/family.html) object. The
  default Gaussian identity family uses the ordinary linear mixed-model
  implementation. Other families are fitted by penalized
  quasi-likelihood (PQL): each outer iteration forms an IRLS working
  response and weights, then fits the resulting weighted Gaussian mixed
  model with the Henderson solver. The initial implementation supports
  one numeric response. Binomial and Poisson working residual
  dispersions are fixed to one. Formula offsets are included on the link
  scale. A supplied symmetric positive definite `W` is combined with the
  IRLS precision after observation filtering, and may be sparse and
  non-diagonal.

- pqlControl:

  A named list controlling non-Gaussian PQL fits. Supported entries are
  `maxit` (outer-iteration limit, default 20) and `tol` (relative
  deviance convergence tolerance, default `1e-5`).

- .pqlInner:

  Internal-use logical indicating that `mmes()` is being called for a
  Gaussian working-model fit within the PQL iteration. Users should not
  set this argument directly.

- .pqlFixedDispersion:

  Internal-use logical indicating that the working residual dispersion
  is fixed, as required for binomial and Poisson PQL models. Users
  should not set this argument directly.

- .pqlWorkingPrecision:

  Internal-use vector containing the current IRLS working precisions
  after observation filtering. Users should not set this argument
  directly.

- .pqlBaseW:

  Internal-use base observation precision matrix supplied by the outer
  PQL fit before combination with the IRLS working precisions. Users
  should not set this argument directly.

- .pqlBaseFactor:

  Internal-use cached factor of `.pqlBaseW`, reused across PQL
  iterations. Users should not set this argument directly.

## Details

**Evaluation environments and observation filtering**

The current interface follows R-style model evaluation. Variables
referenced by `fixed`, `random`, and `rcov` may be supplied as columns
of `data` or resolved from the formula/calling environment. This makes
it possible, for example, to keep relationship matrices, external
grouping variables, or spatial coordinates outside the analysis data
frame.

Missingness is handled through a common observation-selection mechanism.
The response, fixed effects, random effects, and residual covariance
specification are evaluated against the original observation set and
their validity is combined into a single inclusion mask. The same
retained rows are then used for the response vector, fixed-effect design
matrix, random-effect design matrices, optional weight matrix, and
residual covariance indexing. This avoids independently subsetting
different parts of the mixed model.

The fitted object contains `obsInfo`, an observation-level audit table
describing the filtering process. It records the original row index,
validity of the response, fixed, random, and residual components,
whether the observation was retained, and, when applicable, the reason
for exclusion. This can be useful for diagnosing model specifications
involving missing values.

**Estimated marginal means**

When the suggested emmeans package is installed, univariate `mmes` fits
can be supplied directly to
[`emmeans::emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.html).
The resulting means are population-level fixed-effect marginal means
with covariance obtained from the mixed-model equations. Inference uses
asymptotic degrees of freedom (`df = Inf`). Multivariate-response fits
are not currently supported by this interface.

Random and residual formulas are evaluated as R language objects rather
than by splitting their textual representation at plus signs.
Consequently, arithmetic or nested expressions inside model terms are
not interpreted as separate top-level random or residual terms.

The use of this function requires a good understanding of mixed models.
Please review the 'sommer.quick.start' vignette and pay attention to
details like format of your random and fixed variables (e.g. character
and factor variables have different properties when returning BLUEs or
BLUPs).

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

[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md)`(`[`atm`](https://covaruber.github.io/sommer/reference/atm.md)`(x,levels),ism(y))`

can be used to specify heterogeneous variance for the "y" covariate at
specific levels of the covariate "x", e.g.,
*random=~vsm(at(Location,c("A","B")),ism(ID))* fits a variance component
for ID at levels A and B of the covariate Location.

[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md)`(`[`dsm`](https://covaruber.github.io/sommer/reference/dsm.md)`(x),ism(y))`

can be used to specify a diagonal covariance structure for the "y"
covariate for all levels of the covariate "x", e.g.,
*random=~vsm(dsm(Location),ism(ID))* fits a variance component for ID at
all levels of the covariate Location.

[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md)`(`[`usm`](https://covaruber.github.io/sommer/reference/usm.md)`(x),ism(y))`

can be used to specify an unstructured covariance structure for the "y"
covariate for all levels of the covariate "x", e.g.,
*random=~vsm(usm(Location),ism(ID))* fits variance and covariance
components for ID at all levels of the covariate Location.

[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md)`(`[`usm`](https://covaruber.github.io/sommer/reference/usm.md)`(`[`rrm`](https://rdrr.io/pkg/enhancer/man/rrm.html)`(x,y,z,nPC)),ism(y))`

can be used to specify an unstructured covariance structure for the "y"
effect for all levels of the covariate "x", and a response variable "z",
e.g., *random=~vsm(rrm(Location,ID,response, nPC=2),ism(ID))* fits a
reduced-rank factor analytic covariance for ID at 2 principal components
of the covariate Location.

[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md)`(ism(`[`overlay`](https://rdrr.io/pkg/enhancer/man/overlay.html)`(...,rlist=NULL,prefix=NULL)))`

can be used to specify overlay of design matrices between consecutive
random effects specified, e.g., *random=~vsm(ism(overlay(male,female)))*
overlays (overlaps) the incidence matrices for the male and female
random effects to obtain a single variance component for both effects.
The \`rlist\` argument is a list with each element being a numeric value
that multiplies the incidence matrix to be overlayed. See
[`overlay`](https://rdrr.io/pkg/enhancer/man/overlay.html) for
details.Can be combined with vsm().

[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md)`(ism(`[`redmm`](https://rdrr.io/pkg/enhancer/man/redmm.html)`(x,M,nPC)))`

can be used to create a reduced model matrix of an effect (x) assumed to
be a linear function of some feature matrix (M), e.g.,
*random=~vsm(ism(redmm(x,M)))* creates an incidence matrix from a very
large set of features (M) that belong to the levels of x to create a
reduced model matrix. See
[`redmm`](https://rdrr.io/pkg/enhancer/man/redmm.html) for details.Can
be combined with vsm().

[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md)`(`[`leg`](https://rdrr.io/pkg/enhancer/man/leg.html)`(x,n),ism(y))`

can be used to fit a random regression model using a numerical variable
`x` that marks the trayectory for the random effect `y`. The leg
function can be combined with the special functions `dsm`, `usm` `at`
and `csm`. For example *random=~vsm(leg(x,1),ism(y))* or
*random=~vsm(usm(leg(x,1)),ism(y))*.

[`spl2Dc`](https://covaruber.github.io/sommer/reference/spl2Dc.md)`(x.coord, y.coord, at.var, at.levels))`

can be used to fit a 2-dimensional spline (e.g., spatial modeling) using
coordinates `x.coord` and `y.coord` (in numeric class) assuming multiple
variance components. The 2D spline can be fitted at specific levels
using the `at.var` and `at.levels` arguments. For example
*random=~spl2Dc(x.coord=Row,y.coord=Range,at.var=FIELD)*.

**Covariance between random effects**

[`covm`](https://covaruber.github.io/sommer/reference/covm.md)`( `[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md)`(`[`ism`](https://covaruber.github.io/sommer/reference/ism.md)`(ran1)), `[`vsm`](https://covaruber.github.io/sommer/reference/vsm.md)`(`[`ism`](https://covaruber.github.io/sommer/reference/ism.md)`(ran2)) )`

can be used to specify covariance between two different random effects,
e.g., *random=~covm( vsm(ism(x1)), vsm(ism(x2)) )* where two random
effects in their own vsm() structure are encapsulated. Only applies for
simple random effects.

**S3 methods**

S3 methods are available for some parameter extraction such as
[`fitted.mmes`](https://covaruber.github.io/sommer/reference/fitted_mmes.md),
[`residuals.mmes`](https://covaruber.github.io/sommer/reference/residuals_mmes.md),
[`summary.mmes`](https://covaruber.github.io/sommer/reference/summary_mmes.md),
[`randef`](https://covaruber.github.io/sommer/reference/randef.md),
[`coef.mmes`](https://covaruber.github.io/sommer/reference/coef_mmes.md),
[`anova.mmes`](https://covaruber.github.io/sommer/reference/anova_mmes.md),
[`plot.mmes`](https://covaruber.github.io/sommer/reference/plot_mmes.md),
and
[`predict.mmes`](https://covaruber.github.io/sommer/reference/predict_mmes.md)
to obtain adjusted means. In addition, the
[`vpredict`](https://covaruber.github.io/sommer/reference/vpredict.md)
function (replacement of the pin function) can be used to estimate
standard errors for linear combinations of variance components (e.g.,
ratios like h2). The
[`r2`](https://covaruber.github.io/sommer/reference/r2.md) function
calculates reliability.

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

- data:

  the dataset used in model fitting after application of the common
  observation mask.

- obsInfo:

  an observation-level audit table describing the common missing-data
  filter. It includes the original row index, validity indicators for
  response, fixed, random, and residual model components, the final
  inclusion indicator, and the reason for exclusion when applicable.

- Dtable:

  the table to be used for the predict function to help the program
  recognize the factors available.

- llik:

  the vector of log-likelihoods across iterations

- b:

  the vector of fixed effect.

- u:

  the vector of random effect.

- bu:

  the vector of fixed and random effects together.

- Ci:

  the inverse of the coefficient matrix.

- Ci_11:

  the inverse of the coefficient matrix pertaining to the fixed effects.

- theta:

  a list of estimated variance covariance matrices. Each element of the
  list corresponds to the different random and residual components

- covPar:

  the fitted covariance parameters in their descriptor-defined reported
  coordinates.

- covParNative:

  a data frame containing all fixed and estimated covariance parameters
  in the native model-facing scale defined by each CovarianceFactor. For
  example,
  [`dsm()`](https://covaruber.github.io/sommer/reference/dsm.md) entries
  are level-specific variances rather than variance ratios. This is the
  same result returned by
  [`covparams_mmes()`](https://covaruber.github.io/sommer/reference/covparams_mmes.md).

- covParNativeSE:

  the native-scale covariance-parameter table augmented with
  delta-method standard errors and Z ratios. This is the same result
  returned by
  [`covparams_mmes_se()`](https://covaruber.github.io/sommer/reference/covparams_mmes_se.md).

- theta_se:

  inverse of the information matrix.

- InfMat:

  information matrix.

- monitor:

  The values of the variance-covariance components across iterations
  during the REML estimation.

- AIC:

  Akaike information criterion

- BIC:

  Bayesian information criterion

- convergence:

  a TRUE/FALSE statement indicating if the model converged.

- partitions:

  a list where each element contains a matrix indicating where each
  random effect starts and ends.

- partitionsX:

  a list where each element contains a matrix indicating where each
  fixed effect starts and ends.

- percDelta:

  the matrix of percentage change in deltas (see tolParConvNorm
  argument).

- normMonitor:

  the matrix of the three norms calculated (see tolParConvNorm
  argument).

- toBoundary:

  the matrix of variance components that were forced to the boundary
  across iterations.

- Cchol:

  the Cholesky decomposition of the coefficient matrix.

- y:

  the response vector.

- W:

  the column binded matrix W = \[X Z\]

- uList:

  a list containing the BLUPs in data frame format where rows are levels
  of the random effects and column the different factors at which the
  random effect is fitted. This is specially useful for diagonal and
  unstructured models.

- uPevList:

  Prediction error variances for the corresponding entries of `uList`,
  organized by random term and covariance level. These are variances,
  not BLUPs or standard errors. Use `postPEV()` to calculate them after
  fitting and `postVarU()` for the distinct sampling variances of BLUPs.

- randomPrecision:

  Original-level relationship precisions retained by random term for
  post-fit Henderson VarU calculations, avoiding dependence on external
  `Gu` objects. They define each prior covariance together with the fitted
  covariance-coordinate matrix: $G_k=\Sigma_k\otimes A_k$.
  `postVarU()` adds model-effect outputs, whereas `predict(PEV=TRUE,
  VarU=TRUE)` projects covariances through the requested `D`/`Dtable`
  linear combinations.

- args:

  the fixed, random and residual formulas from the mmes model.

- constraints:

  The vector of constraints.

## References

Covarrubias-Pazaran G. Genome assisted prediction of quantitative traits
using the R package sommer. PLoS ONE 2016, 11(6):
doi:10.1371/journal.pone.0156744

Jensen, J., Mantysaari, E. A., Madsen, P., and Thompson, R. (1997).
Residual maximum likelihood estimation of (co) variance components in
multivariate mixed linear models using average information. Journal of
the Indian Society of Agricultural Statistics, 49, 215-236.

Sanderson, C., & Curtin, R. (2025). Armadillo: An Efficient Framework
for Numerical Linear Algebra. arXiv preprint arXiv:2502.03000.

Gilmour et al. 1995. Average Information REML: An efficient algorithm
for variance parameter estimation in linear mixed models. Biometrics
51(4):1440-1450.

## Author

Coded by Giovanny Covarrubias-Pazaran with contributions of Christelle
Fernandez Camacho, Johan Aparicio-Arce, and Claudio Flavio 5.1 jaja
(laughing in Spanish).

## Examples

``` r
####=========================================####
#### For CRAN time limitations most lines in the
#### examples are silenced with one '#' mark,
#### remove them and run the examples
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
ans1 <- mmes(Yield~Env,
             random= ~ Name + Env:Name,
             rcov= ~ units,
             data=DT)
#> Solver selected: ldlt
#> OpenMP available: up to 8 threads.
#> OpenMP active: parallel SLQ log-determinant probes (8 probes).
#> iteration    LogLik     wall    cpu(sec)   restrained   EM weight      pivot
#>     1      -31.6987   21:18:5      0           0      1      3.08144
#>     2      -28.0631   21:18:5      0           0      0.813615      5.42039
#>     3      -26.3167   21:18:5      0           0      0.661969      3.92118
#>     4      -26.2862   21:18:5      0           0      0.538588      4.00163
#>     5      -26.2851   21:18:5      0           0      0.438203      3.9582
#>     6      -26.2851   21:18:5      0           0      0.356529      3.95514
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

# \donttest{

####===========================================####
#### Univariate heterogeneous variance models  ####
####===========================================####
DT=DT[with(DT, order(Env)), ]
## Compound simmetry (CS) + Diagonal (DIAG) model
ans2 <- mmes(Yield~Env,
             random= ~Name + vsm(dsm(Env),ism(Name)),
             rcov= ~ vsm(dsm(Env),ism(units)),
             data=DT)
#> Solver selected: ldlt
#> OpenMP available: up to 8 threads.
#> OpenMP active: parallel SLQ log-determinant probes (8 probes).
#>   Working-coordinate trust scaling applied: alpha=0.650987
#> iteration    LogLik     wall    cpu(sec)   restrained   EM weight      pivot
#>     1      -31.6987   21:18:5      0           0      1      3.08144
#>     2      -29.0398   21:18:5      0           0      0.813615      4.45068
#>     3      -24.0636   21:18:5      0           0      0.661969      2.4284
#>     4      -22.7208   21:18:5      0           0      0.538588      3.22302
#>     5      -21.8476   21:18:5      0           0      0.438203      2.6077
#>     6      -21.7347   21:18:5      0           0      0.356529      2.54008
#>     7      -21.667   21:18:5      0           0      0.290077      2.4098
#>     8      -21.616   21:18:5      0           0      0.236011      2.2611
#>     9      -21.5897   21:18:5      0           0      0.192022      2.15979
#>     10      -21.5774   21:18:5      0           0      0.156232      2.09536
#>     11      -21.5723   21:18:5      0           0      0.127113      2.05671
#>     12      -21.5703   21:18:5      0           0      0.103421      2.03495
#>     13      -21.5697   21:18:5      0           0      0.0841447      2.02349
#>     14      -21.5696   21:18:5      0           0      0.0684614      2.01786
#>     15      -21.5695   21:18:5      0           0      0.0557012      2.01529
summary(ans2)
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

####===========================================####
####  Univariate unstructured variance models  ####
####===========================================####

ans3 <- mmes(Yield~Env,
             random=~ vsm(usm(Env),ism(Name)),
             rcov=~vsm(dsm(Env),ism(units)),
             data=DT)
#> Solver selected: ldlt
#> OpenMP available: up to 8 threads.
#> OpenMP active: parallel SLQ log-determinant probes (8 probes).
#>   Working-coordinate trust scaling applied: alpha=0.872325
#> iteration    LogLik     wall    cpu(sec)   restrained   EM weight      pivot
#>     1      -32.6293   21:18:5      0           0      1      3.11257
#>     2      -25.0808   21:18:5      0           0      0.813615      2.42213
#>     3      -20.5724   21:18:5      0           0      0.661969      2.11889
#>     4      -19.1091   21:18:5      0           0      0.538588      2.40562
#>     5      -18.394   21:18:5      0           0      0.438203      2.31081
#>     6      -18.1453   21:18:5      0           0      0.356529      2.35239
#>     7      -17.9046   21:18:5      0           0      0.290077      2.37242
#>     8      -17.7402   21:18:5      0           0      0.236011      2.37925
#>     9      -17.6754   21:18:5      0           0      0.192022      2.38326
#>     10      -17.6517   21:18:5      0           0      0.156232      2.38489
#>     11      -17.6433   21:18:5      0           0      0.127113      2.38555
#>     12      -17.6405   21:18:5      0           0      0.103421      2.38578
#>     13      -17.6397   21:18:5      0           0      0.0841447      2.38588
#>     14      -17.6395   21:18:5      0           0      0.0684614      2.38595
#>     15      -17.6394   21:18:5      0           0      0.0557012      2.386
summary(ans3)
#> ============================================================
#>          Multivariate Linear Mixed Model fit by  REML         
#> **********************  sommer 4.4  ********************** 
#> ============================================================
#>         logLik      AIC      BIC Method Converge
#> Value -17.6394 41.27881 50.93988     AI     TRUE
#> ============================================================
#> Variance-Covariance components:
#>                        term factor                   parameter estimate
#> 1  vsm(usm(Env), ism(Name))     us           variance[CA.2011]  15.6594
#> 2  vsm(usm(Env), ism(Name))     us covariance[CA.2012,CA.2011]   6.1105
#> 3  vsm(usm(Env), ism(Name))     us covariance[CA.2013,CA.2011]   6.3837
#> 4  vsm(usm(Env), ism(Name))     us           variance[CA.2012]   4.5320
#> 5  vsm(usm(Env), ism(Name))     us covariance[CA.2013,CA.2012]   0.3913
#> 6  vsm(usm(Env), ism(Name))     us           variance[CA.2013]   8.5980
#> 7 vsm(dsm(Env), ism(units))   diag           variance[CA.2011]   4.9770
#> 8 vsm(dsm(Env), ism(units))   diag           variance[CA.2012]   5.6707
#> 9 vsm(dsm(Env), ism(units))   diag           variance[CA.2013]   2.5566
#>   StdError Zratio
#> 1   3.4580 4.5285
#> 2   1.6633 3.6736
#> 3   1.9750 3.2323
#> 4   1.3093 3.4615
#> 5   1.0846 0.3608
#> 6   1.7656 4.8698
#> 7   0.9077 5.4828
#> 8   0.9161 6.1900
#> 9   0.4516 5.6615
#> ============================================================
#> Fixed effects:
#>            Estimate Std.Error t.value
#> Intercept    16.331        NA      NA
#> EnvCA.2012   -5.696        NA      NA
#> EnvCA.2013   -6.271        NA      NA
#> ============================================================
#> Use the '$' sign to access results and parameters


# }
```
