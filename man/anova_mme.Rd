\name{anova.mme}
\alias{anova.mme}
\title{anova form a GLMM fitted with mme}
\description{\code{anova} method for class \code{"mme"}.
}
\usage{
\method{anova}{mme}(object, object2=NULL, \dots)
}
\arguments{
 \item{object}{an object of class \code{"mme"}}
 \item{object2}{an object of class \code{"mme"}, if NULL the program will provide regular sum of squares results.}
 \item{\dots}{Further arguments to be passed}
}
\value{
  vector of anova
 }
\author{Giovanny Covarrubias}
\seealso{\code{\link{anova}}, \code{\link{mme}}}
\keyword{models}

