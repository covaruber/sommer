\name{predict.mmes}
\alias{predict.mmes}
\alias{predict.mmer}
\title{Predict form of a LMM fitted with mmes}
\description{\code{predict} method for class \code{"mmes"}.
}
\usage{
\method{predict}{mmes}(object, Dtable=NULL, D, levels=NULL, sed=FALSE,
        pairwise=FALSE, adjust="none", df=Inf, \dots)
}
\arguments{
 \item{object}{a mixed model of class \code{"mmes"}}
 \item{Dtable}{a table specifying the terms to be included or averaged.
 
  An "include" term means that the model matrices for that fixed or random effect is filled with 1's for the positions where column names and row names match.
  
  An "include and average" term means that the included cells of that effect are averaged: each prediction row is divided by the number of cells of the term present with it (for fixed terms this count includes the reference cells absorbed in the intercept).
  
  An "average" term alone means that all rows for such fixed or random effect will be filled with 1/#levels in the effect.
  
  If a term is not considered "include" or "average" is then totally ignored in the BLUP and SE calculation. 
  
  The default rule to invoke when the user doesn't provide the Dtable is to include and average all terms that match the argument D.

  The \code{levels} column (a list, \code{NULL} meaning all levels) restricts
  each term to the given levels; see the \code{levels} argument.
 
 }
 \item{D}{a character string specifying the variable used to extract levels for the rows of the D matrix and its construction. Alternatively, the D matrix (of class dgCMatrix) specifying the matrix to be used for the predictions directly.}
 \item{levels}{an optional named list whose names are values of \code{Dtable$term}.
  Each element is stored in the \code{levels} column of the returned Dtable,
  so the Dtable remains the complete description of how D was built. Levels
  are written as data levels for fixed factors (\code{"a:b"} for
  interactions, e.g. \code{"Victory:0.2"}), as coefficient labels for random
  terms, and as numeric values for numeric covariate terms (the value at which
  the covariate is predicted). For an included term matching the classify
  variable(s) the levels select the prediction rows, in the given order, and
  may include levels without records (e.g. levels only present in \code{Gu}).
  For other included or averaged terms, only the listed levels enter the
  inclusion or the average.}
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
 \item{\dots}{Further arguments to be passed.}
}
\details{
  This function allows to produce predictions specifying those variables that define the margins of the hypertable to be predicted (argument D). Predictions are obtained for each combination of values of the specified variables that is present in the data set used to fit the model. See vignettes for more details. 
  
  Standard errors are exact for any \code{computeCi} setting (0, 1, or 2) and do not require \code{object$Ci}: \code{Var(D \%*\% bu) = D C.inv() D.t()} is instead obtained by solving \code{C x = d} for each row \code{d} of \code{D} against the stored coefficient matrix \code{object$C}, which \code{mmes()} always returns regardless of \code{computeCi}. This is cheaper than forming the complete inverse (\code{computeCi=2}) and, unlike the Takahashi selected-inverse subset (\code{computeCi=1}), is exact for arbitrary linear combinations, not just diagonal PEVs.
  
  For predicted values the pertinent design matrices X and Z together with BLUEs (b) and BLUPs (u) are multiplied and added together.
  
  predicted.value equal Xb + Zu.1 + ... + Zu.n
  
  For computing standard errors for predictions the parts of the coefficient matrix:
  
  C11 equal (X.t() V.inv() X).inv()
  
  C12 equal 0 - [(X.t() V.inv() X).inv() X.t() V.inv() G Z] 
  
  C22 equal PEV equal G - [Z.t() G[V.inv() - (V.inv() X X.t() V.inv() X V.inv() X)]G Z.t()]
  
  In practive C equals ( W.t() V.inv() W ).inv()
  
  when both fixed and random effects are present in the inclusion set. If only fixed and random effects are included, only the respective terms from the SE for fixed or random effects are calculated.
  
 }
\value{

  \item{pvals}{the table of predictions according to the specified arguments.}
  
  \item{vcov}{the variance covariance for the predictions.}
  
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

}
