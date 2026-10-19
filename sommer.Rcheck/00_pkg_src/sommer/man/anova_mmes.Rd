\name{anova.mmes}
\alias{anova.mmes}
\alias{anova.mmer}
\title{anova form a GLMM fitted with mmes}
\description{\code{anova} method for class \code{"mmes"}.
}
\usage{
\method{anova}{mmes}(object, object2=NULL, \dots)
}
\arguments{
 \item{object}{an object of class \code{"mmes"}}
 \item{object2}{an object of class \code{"mmes"}. If supplied, a likelihood ratio test between the two models is returned (not available for PQL fits). If \code{NULL}, Wald tests for the fixed effects of \code{object} are returned by \code{\link{wald_mmes}}.}
 \item{\dots}{Further arguments passed to \code{\link{wald_mmes}} (e.g. \code{ssType}, \code{denDF}) when \code{object2} is \code{NULL}.}
}
\value{
  vector of anova
 }
\author{Giovanny Covarrubias}
\seealso{\code{\link{anova}}, \code{\link{mmes}}, \code{\link{wald_mmes}}}
\keyword{models}

