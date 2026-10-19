\name{predict.mmes}
\alias{predict.mmes}
\alias{predict.mmer}
\title{Predict form of a LMM fitted with mmes}
\description{\code{predict} method for class \code{"mmes"}.
}
\usage{
\method{predict}{mmes}(object, Dtable=NULL, D, levels=NULL, sed=FALSE,
        pairwise=FALSE, adjust="none", df=Inf, PEV=TRUE, VarU=TRUE, \dots)
}
\arguments{
 \item{object}{a mixed model of class \code{"mmes"}}
 \item{Dtable}{a table specifying the terms to be included or averaged.
 
  An "include" term uses the cells matching each prediction row. For
  categorical random effects these are mapped within the fitted term, including
  fitted cells without records. Numeric fixed terms retain their evaluated
  covariate values rather than being replaced by indicators.
  
  An "include and average" term means that the included cells of that effect are averaged: each prediction row is divided by the number of cells of the term present with it (for fixed terms this count includes the reference cells absorbed in the intercept).
  
  An "average" term alone averages over the selected factor levels. For
  categorical terms the weights are 1/#levels; numeric fixed terms are
  evaluated at their covariate settings before averaging.
  
  If a term is not considered "include" or "average" is then totally ignored in the BLUP and SE calculation. 
  
  The default rule to invoke when the user doesn't provide the Dtable is to include and average all terms that match the argument D.

  The \code{levels} column is a list-column. \code{NULL} means all factor
  levels and default numeric settings. Entries restrict levels or specify
  covariate settings; see the \code{levels} argument. The returned Dtable also
  records resolved settings for numeric terms and can be reused.
 
 }
 \item{D}{a character string specifying the variable used to extract levels for the rows of the D matrix and its construction. Alternatively, the D matrix (of class dgCMatrix) specifying the matrix to be used for the predictions directly.}
 \item{levels}{an optional named list whose names are values of \code{Dtable$term}.
  Each element is stored in the \code{levels} column of the returned Dtable,
  so the Dtable remains the complete description of how D was built. Levels
  are written as data levels for fixed factors (\code{"a:b"} for
  interactions, e.g. \code{"Victory:0.2"}), as cell labels for categorical random
  terms (e.g. \code{"CA.2011:genotype1"} for an Env by Name term), and as numeric
  settings for covariate terms. A term containing one numeric variable accepts
  a numeric scalar, e.g. \code{levels=list("Env:cov"=5)}. Numeric or mixed
  terms also accept a named per-variable list, e.g.
  \code{levels=list("Env:cov"=list(Env=c("CA.2011", "CA.2012"), cov=5))}.
  Each numeric setting must be one finite value. Use the exact term spelling
  in \code{object$Dtable$term}. For an included categorical term matching the classify
  variable(s) the levels select the prediction rows, in the given order, and
  may include levels without records (e.g. levels only present in \code{Gu}).
  For other included or averaged terms, only the listed levels enter the
  inclusion or the average. Restrictions are term-specific: selecting levels
  of \code{V} does not automatically restrict \code{V:N}.}
 \item{sed}{if \code{TRUE}, return the matrix of standard errors of differences
  \code{sed} and its summary \code{avsed} (mean, minimum and maximum).}
 \item{pairwise}{\code{TRUE} for all pairwise differences between predictions,
  or a single prediction level to compare every other prediction against it.}
 \item{adjust}{multiplicity adjustment passed to \code{\link[stats]{p.adjust}}.}
 \item{df}{degrees of freedom for pairwise tests: \code{Inf} (normal-based,
  default), a positive number, \code{"residual"} for n - rank(X), or
  \code{"satterthwaite"}/\code{"kr"} (see \code{\link{wald_mmes}}). The latter
  two apply to differences involving only fixed effects; differences involving
  random effects keep normal-based tests (\code{df=Inf}). With \code{"kr"} the
  SED of those differences uses the Kenward-Roger adjusted covariance.}
 \item{PEV}{Logical, default \code{TRUE}. Return the full target prediction-error
 covariance as \code{PEV} and its diagonal as \code{pvals$pev}. When
 \code{FALSE}, only these additional outputs are omitted; \code{vcov} and
 \code{std.error} remain available for compatibility.}
 \item{VarU}{Logical, default \code{TRUE}. Return the sampling covariance of
 the random contribution as \code{VarU}, its diagonal as \code{pvals$var.u},
 and the sampling covariance of the entire predictor as \code{sampling.vcov}.
 Set \code{FALSE} to skip the additional relationship solves and outputs.
 These calculations require a Gaussian Henderson fit. For PQL predictions
 use \code{VarU=FALSE}; their existing working-model uncertainty remains an
 approximation.}
 \item{\dots}{Further arguments to be passed.}
}
\details{
  This function allows to produce predictions specifying those variables that define the margins of the hypertable to be predicted (argument D). Predictions are obtained for each combination of values of the specified variables that is present in the data set used to fit the model. See vignettes for more details. 

  For categorical structured random terms, the fitted covariance-coordinate
  and main-effect levels define the genetic prediction weights, even when a
  cell has no observations. The \code{levels} argument can additionally request
  fitted levels present only in \code{Gu}. Unsupported random designs that
  cannot be mapped for all prediction rows produce a warning; an explicit
  \code{D} matrix can be used for these designs.

  Numeric fixed terms, including factor-by-covariate interactions, are
  evaluated using the fitted fixed formula and contrasts. Unless overridden,
  each numeric variable is set to its mean in the retained fitting data.
  Factor levels are equally weighted when averaged. Thus an averaged
  \code{Env:cov} term with one slope per environment uses
  \code{mean(cov)/nlevels(Env)} for each slope, not \code{1/nlevels(Env)}.
  
  Write \eqn{D=[L_b,L_u]}. The target is \eqn{t=L_b\beta+L_u u} and its
  predictor is \eqn{\hat t=L_b\hat\beta+L_u\hat u}. The existing
  prediction \code{PEV} and \code{VarU} matrices are for the linear
  combinations represented by the final returned \code{D}, after applying
  \code{Dtable} and \code{levels}, not copies of model-level effect variances.
  Each covariance matrix has one row and column per requested prediction.
  The existing
  \code{vcov}, optional \code{PEV}, and squared \code{std.error} refer to
  \eqn{\operatorname{Var}(t-\hat t)}, not generally to
  \eqn{\operatorname{Var}(\hat t)}. They are computed from
  \deqn{\operatorname{Var}(t-\hat t)=D C_0^{-1}D',}
  where \eqn{C_0} is the original-scale Henderson coefficient matrix defined
  in \code{\link{postVarU}}. The stored system \code{C} is solved against
  contrast transposes and the result multiplied by \code{Cscale}. Arbitrary
  contrasts work regardless of \code{computeCi}; no stored \code{Ci} is used.

  With \eqn{Q=(C_0^{-1})_{uu}} and \eqn{B=(C_0^{-1})_{bb}}, the new outputs are
  \deqn{\mathrm{VarU}=L_u(G-Q)L_u',}
  \deqn{\mathrm{sampling.vcov}=L_b B L_b'+L_u(G-Q)L_u'.}
  The BLUE and BLUP have zero sampling cross-covariance, whereas the
  fixed/random blocks of the inverse Henderson matrix describe error
  covariance and generally are not zero. In particular, subtracting the full
  target PEV from \eqn{L_uGL_u'} is incorrect if fixed effects are included.
  A fixed-only prediction has zero \code{VarU}, but generally nonzero
  \code{sampling.vcov} and PEV. Random-only predictions satisfy
  \eqn{L_uGL_u'=\mathrm{PEV}+\mathrm{VarU}}.

  The implementation batches full, random-only and fixed-only contrasts into
  one Henderson factorization. Prior products use the term structure
  \eqn{G_k=\Sigma_k\otimes A_k} and sparse relationship-precision solves.
  Neither the observation covariance \eqn{V} nor its inverse is constructed.
  Additional matrices require quadratic storage in the number of predictions.
  Rotation and factor-score mappings operate in public effect coordinates.
  New fits retain relationship precisions; older ordinary fits attempt model
  input reconstruction as described in \code{\link{postVarU}}.

  All covariance outputs use the fitted variance components as fixed,
  omitting their estimation uncertainty. They are exact for known components
  and plug-in approximations for estimated components, for both ML and REML.
  SEDs, pairwise tests and \code{std.error} retain their PEV-based meaning;
  enabling VarU does not change those calculations. Direct and matrix-free
  \code{solveOnly} fits lack the stored coefficient system required here.
  
 }
\value{

  \item{pvals}{the table of predictions according to the specified arguments.}
  
  \item{vcov}{The full target prediction-error covariance, always returned.}
  \item{PEV}{When requested, identical to \code{vcov}; its diagonal is
  \code{pvals$pev}.}
  \item{VarU}{When requested, sampling covariance of the predicted random
  contribution only; its diagonal is \code{pvals$var.u}.}
  \item{sampling.vcov}{When \code{VarU=TRUE}, sampling covariance of the entire
  fixed-plus-random predictor. This is not the covariance used for standard errors.}
  
  \item{D}{the model matrix for predictions as defined in Welham et al.(2004).}
  
  \item{Dtable}{the table specifying the terms to include and terms to be averaged, including the \code{levels} used.}

  \item{sed, avsed}{when \code{sed=TRUE}: the standard errors of differences, sqrt(v_ii + v_jj - 2 v_ij) from \code{vcov}, and their mean, minimum and maximum.}

  \item{pairwise}{when requested: a data frame with \code{level1}, \code{level2}, \code{difference} (level1 - level2), \code{SED}, \code{statistic}, \code{df}, \code{p.value} and \code{p.adjusted}.}
  
 }
 \references{
  Welham, S., Cullis, B., Gogel, B., Gilmour, A., and Thompson, R. (2004). Prediction in linear mixed models. Australian and New Zealand Journal of Statistics, 46, 325 - 347.
 }
\author{Giovanny Covarrubias-Pazaran}
\seealso{\code{\link{predict}}, \code{\link{mmes}}}
\keyword{models}
\examples{

data(DT_yatesoats, package="enhancer")
DT <- DT_yatesoats
m3 <- mmes(fixed=Y ~ V + N + V:N ,
           random = ~ B + B:MP,
           rcov=~units,
           data = DT)

m3 <- postPEV(m3, mode = 2)
 
#############################
## predict means for nitrogen
#############################
Dt <- m3$Dtable; Dt
# first fixed effect just average
Dt[1,"average"] = TRUE
# second fixed effect include
Dt[2,"include"] = TRUE
# third fixed effect include and average
Dt[3,"include"] = TRUE
Dt[3,"average"] = TRUE
Dt

pp=predict(object=m3, Dtable=Dt, D="N")
pp$pvals

#############################
## predict means for variety
#############################

Dt <- m3$Dtable; Dt
# first fixed effect include
Dt[1,"include"] = TRUE
# second fixed effect just average
Dt[2,"average"] = TRUE
# third fixed effect include and average
Dt[3,"include"] = TRUE
Dt[3,"average"] = TRUE
Dt

pp=predict(object=m3, Dtable=Dt, D="V")
pp$pvals

#############################
## predict means for nitrogen:variety
#############################
# prediction matrix D based on (equivalent to classify in asreml)
Dt <- m3$Dtable; Dt
# first fixed effect include and average
Dt[1,"include"] = TRUE
Dt[1,"average"] = TRUE
# second fixed effect include and average
Dt[2,"include"] = TRUE
Dt[2,"average"] = TRUE
# third fixed effect include and average
Dt[3,"include"] = TRUE
Dt[3,"average"] = TRUE
Dt

pp=predict(object=m3, Dtable=Dt, D="N:V")
pp$pvals

#############################
## select levels through the argument or the Dtable
#############################
Dt <- m3$Dtable
Dt$average[Dt$term \%in\% c("1", "V", "V:N")] <- TRUE
Dt$include[Dt$term \%in\% c("N", "V:N")] <- TRUE

p <- predict(m3, D="N", Dtable=Dt,
             levels=list(N=c("0.6", "0.2")))
p$pvals
p$Dtable$levels[[match("N", p$Dtable$term)]]

Dt$levels[[match("N", Dt$term)]] <- c("0.6", "0.2")
p2 <- predict(m3, D="N", Dtable=Dt)
all.equal(p$D, p2$D)

## Covariances of the requested linear combinations
p <- predict(m3, D="N", Dtable=Dt, PEV=TRUE, VarU=TRUE)
p$PEV
p$VarU
p$sampling.vcov
randomColumns <- seq.int(nrow(m3$b)+1L, nrow(m3$bu))
randomD <- p$D[,randomColumns,drop=FALSE]
full <- postVarU(m3, mode=2)
all.equal(unname(p$VarU),
          unname(as.matrix(randomD \%*\% full$VarU \%*\% t(randomD))))

}
