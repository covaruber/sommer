\name{residuals.mmec}
\alias{residuals.mmec}
\title{Residuals form a GLMM fitted with mmec}
\description{\code{residuals} method for class \code{"mmec"}.
}
\usage{
\method{residuals}{mmec}(object, \dots)
}
\arguments{
 \item{object}{an object of class \code{"mmec"}}
 \item{\dots}{Further arguments to be passed}
}
\value{
  vector of residuals of the form e = y - Xb - Zu, the so called conditional residuals.
 }
\author{Giovanny Covarrubias}
\seealso{\code{\link{residuals}}, \code{\link{mmec}}}
\keyword{models}

