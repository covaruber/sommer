\name{residuals.mmer}
\alias{residuals.mmer}
\title{Residuals form a GLMM fitted with mmer}
\description{\code{residuals} method for class \code{"mmer"}.
}
\usage{
\method{residuals}{mmer}(object, \dots)
}
\arguments{
 \item{object}{an object of class \code{"mmer"}}
 \item{\dots}{Further arguments to be passed}
}
\value{
  vector of residuals of the form e = y - Xb - Zu, the so called conditional residuals.
 }
\author{Giovanny Covarrubias}
\seealso{\code{\link{residuals}}, \code{\link{mmer}}}
\keyword{models}

