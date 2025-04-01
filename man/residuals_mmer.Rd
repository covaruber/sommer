\name{residuals.mme}
\alias{residuals.mme}
\title{Residuals form a GLMM fitted with mme}
\description{\code{residuals} method for class \code{"mme"}.
}
\usage{
\method{residuals}{mme}(object, \dots)
}
\arguments{
 \item{object}{an object of class \code{"mme"}}
 \item{\dots}{Further arguments to be passed}
}
\value{
  vector of residuals of the form e = y - Xb - Zu, the so called conditional residuals.
 }
\author{Giovanny Covarrubias}
\seealso{\code{\link{residuals}}, \code{\link{mme}}}
\keyword{models}

