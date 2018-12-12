\name{anova.mmer}
\alias{anova.mmer}
\title{anova form a GLMM fitted with mmer}
\description{\code{anova} method for class \code{"mmer"}.
}
\usage{
\method{anova}{mmer}(object, object2=NULL, type=1, \dots)
}
\arguments{
 \item{object}{an object of class \code{"mmer"}}
 \item{object2}{an object of class \code{"mmer"}, if NULL the program will provide regular sum of squares results.}
  \item{type}{anova type, I or II}
 \item{\dots}{Further arguments to be passed}
}
\value{
  vector of anova
 }
\author{Giovanny Covarrubias \email{covarrubiasp@wisc.edu}}
\seealso{\code{\link{anova}}, \code{\link{mmer}}}
\keyword{models}

