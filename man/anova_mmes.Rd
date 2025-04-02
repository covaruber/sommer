\name{anova.mmes}
\alias{anova.mmes}
\title{anova form a GLMM fitted with mmes}
\description{\code{anova} method for class \code{"mmes"}.
}
\usage{
\method{anova}{mmes}(object, object2=NULL, \dots)
}
\arguments{
 \item{object}{an object of class \code{"mmes"}}
 \item{object2}{an object of class \code{"mmes"}, if NULL the program will provide regular sum of squares results.}
 \item{\dots}{Further arguments to be passed}
}
\value{
  vector of anova
 }
\author{Giovanny Covarrubias}
\seealso{\code{\link{anova}}, \code{\link{mmes}}}
\keyword{models}

