\name{residuals.mmes}
\alias{residuals.mmes}
\alias{residuals.mmer}
\title{Residuals form a GLMM fitted with mmes}
\description{\code{residuals} method for class \code{"mmes"}.
}
\usage{
\method{residuals}{mmes}(object,
                          type=c("response", "deviance", "pearson", "working",
                            "quantile"), \dots)
}
\arguments{
 \item{object}{an object of class \code{"mmes"}}
 \item{type}{For a non-Gaussian PQL fit, return response residuals (default),
 signed deviance residuals, Pearson residuals, or final IRLS working residuals.
 Deviance and Pearson residuals use the prior weights (binomial trials, see
 \code{\link{binm}}). Native Laplace fits also support randomized quantile
 residuals. Discrete responses use randomized probability integral
 transforms; Gaussian and Gamma responses use their fitted continuous CDFs.
 These residuals are conditional on fitted model parameters and latent modes.}
 \item{\dots}{Further arguments to be passed}
}
\value{
  For Gaussian fits, residuals of the form e = y - Xb - Zu. For non-Gaussian
  PQL fits, residuals on the requested scale. Native Laplace quantile residuals
  are returned on the standard-normal scale.
 }
\author{Giovanny Covarrubias}
\seealso{\code{\link{residuals}}, \code{\link{mmes}}}
\keyword{models}
