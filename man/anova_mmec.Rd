\name{anova.mmec}
\alias{anova.mmec}
\title{anova form a GLMM fitted with mmec}
\description{\code{anova} method for class \code{"mmec"}.
}
\usage{
\method{anova}{mmec}(object, object2=NULL, \dots)
}
\arguments{
 \item{object}{an object of class \code{"mmec"}}
 \item{object2}{an object of class \code{"mmec"}, if NULL the program will provide regular sum of squares results.}
 \item{\dots}{Further arguments to be passed}
}
\value{
  vector of anova
 }
\author{Giovanny Covarrubias}
\seealso{\code{\link{anova}}, \code{\link{mmec}}}
\keyword{models}

