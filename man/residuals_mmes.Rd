\name{residuals.mmes}
\alias{residuals.mmes}
\alias{residuals.mmer}
\title{Residuals form a GLMM fitted with mmes}
\description{\code{residuals} method for class \code{"mmes"}.
}
\usage{
\method{residuals}{mmes}(object,
                          type=c("response", "deviance", "working"), \dots)
}
\arguments{
 \item{object}{an object of class \code{"mmes"}}
 \item{type}{For a non-Gaussian PQL fit, return response residuals (default),
 signed deviance residuals, or final IRLS working residuals.}
 \item{\dots}{Further arguments to be passed}
}
\value{
  For Gaussian fits, residuals of the form e = y - Xb - Zu. For non-Gaussian
  PQL fits, residuals on the requested scale.
 }
\author{Giovanny Covarrubias}
\seealso{\code{\link{residuals}}, \code{\link{mmes}}}
\keyword{models}

