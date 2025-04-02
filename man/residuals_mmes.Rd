\name{residuals.mmes}
\alias{residuals.mmes}
\title{Residuals form a GLMM fitted with mmes}
\description{\code{residuals} method for class \code{"mmes"}.
}
\usage{
\method{residuals}{mmes}(object, \dots)
}
\arguments{
 \item{object}{an object of class \code{"mmes"}}
 \item{\dots}{Further arguments to be passed}
}
\value{
  vector of residuals of the form e = y - Xb - Zu, the so called conditional residuals.
 }
\author{Giovanny Covarrubias}
\seealso{\code{\link{residuals}}, \code{\link{mmes}}}
\keyword{models}

