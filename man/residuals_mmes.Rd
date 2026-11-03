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
 These residuals are conditional on fitted model parameters and latent modes.
 Nominal multinomial fits support response, Pearson and unsigned deviance
 residuals, but not working or quantile residuals.}
 \item{\dots}{Further arguments to be passed}
}
\value{
  For Gaussian fits, residuals of the form e = y - Xb - Zu. For non-Gaussian
  PQL fits, residuals on the requested scale. Native Laplace quantile residuals
  are returned on the standard-normal scale.
  For \code{\link{ordinal}} fits, response residuals subtract the expected
  numeric category code; Pearson residuals divide by its conditional standard
  deviation and working residuals by its predictor derivative. These three
  types depend on category coding. Deviance residuals have squared magnitude
  \eqn{-2\log P(Y_i=y_i\mid\hat u)}, with positive sign when the observed
  code is at least the expected code. Quantile residuals randomize over the
  observed category's cumulative-probability interval. Undefined zero
  variance or zero-derivative residuals raise an error.
  For \code{\link{multinomial}} fits, response and Pearson residuals are
  observation-by-category matrices with the original category labels.
  Response residuals are \eqn{1(Y_i=k)-p_{ik}} and Pearson residuals divide
  each cell by \eqn{\sqrt{p_{ik}(1-p_{ik})}}. Category residuals within a row
  are dependent; do not treat them as independent normal observations.
  Deviance residuals are a non-negative vector with magnitude
  \eqn{\sqrt{-2\log p_{i,y_i}}}. No sign is assigned because nominal
  categories have no natural order. Squared magnitudes sum to the fitted
  conditional response deviance, not minus twice the marginal Laplace
  likelihood. Log probabilities are evaluated from logits without clipping,
  so large finite tail contributions remain visible even when probabilities
  underflow. Pearson residuals at probabilities zero or one raise an error.
  Working and quantile residuals remain unsupported.
  With categorical \code{\link{familym}}, residuals are a trait-named list
  applying these definitions independently to each trait's retained rows.
  Scalar and ordinal components are vectors; nominal response/Pearson
  components are labelled matrices. A requested residual type unsupported
  by any declared trait raises an error rather than assigning an arbitrary
  nominal ordering.
  These are descriptive, in-sample model-checking quantities conditional on
  fitted parameters and random modes, not calibrated outlier tests or
  goodness-of-fit p-values. Large deviance magnitudes identify outcomes the
  fitted model considers unlikely; compare category-wise residuals against
  covariates or groups to investigate patterns. No uncertainty calculation
  or refitting is performed by \code{residuals()}.
 }
\examples{
\dontrun{
residuals(fit, type="deviance")  # unsigned nominal lack-of-fit magnitudes
residuals(fit, type="pearson")   # labelled category-wise residual matrix
}
}
\author{Giovanny Covarrubias}
\seealso{\code{\link{residuals}}, \code{\link{mmes}}}
\keyword{models}
