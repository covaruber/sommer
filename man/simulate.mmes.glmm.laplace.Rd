\name{simulate.mmes.glmm.laplace}
\alias{simulate.mmes.glmm.laplace}
\title{Simulate responses from a fitted native Laplace GLMM}
\description{
Generate response replicates from the fitted native family distributions
conditional on the fitted fixed effects and conditional random-effect modes.
For population simulations, new random effects are drawn from their fitted
population covariance using sparse relationship-precision factors.
}
\usage{
\method{simulate}{mmes.glmm.laplace}(object, nsim=1, seed=NULL,
  re.form=c("conditional", "population"), maxMemoryMB=NULL, \dots)
}
\arguments{
  \item{object}{A fitted native Laplace GLMM.}
  \item{nsim}{Number of response replicates to generate; a positive integer.}
  \item{seed}{Optional integer seed. When supplied, the caller's existing
  random-number state is restored after simulation.}
  \item{re.form}{\code{"conditional"} simulates responses using the fitted
  conditional random-effect modes. \code{"population"} first draws new
  Gaussian random effects from the fitted covariance model, then draws
  responses.}
  \item{maxMemoryMB}{Maximum estimated workspace in megabytes for population
  random-effect draws and the returned simulation matrix. The default uses
  \code{execution$denseMemoryMB}; sparse-factor fill is data-dependent.}
  \item{\dots}{Unused; additional arguments are rejected.}
}
\value{
A data frame with one column per replicate and one row per fitted observation.
Grouped-binomial responses are returned as proportions, matching the response
scale used by the fitted object.
}
\details{
Binomial, Poisson, NB2, Gaussian and Gamma draws use the fitted family
parameters, including estimated NB2 shape and continuous-family dispersion.
For long-format mixed-family fits, each row is simulated from its declared
family. Conditional simulations hold the fitted random effects fixed.
Population simulations draw random effects jointly by covariance term using
sparse relationship-precision Cholesky solves and dense covariance-shape
factors, subject to the approximate workspace limit \code{maxMemoryMB}. These
simulations are not a parametric-bootstrap calibration procedure and do not
include parameter-estimation uncertainty.
}
\seealso{\code{\link{mmes}}, \code{\link{residuals.mmes}}}
\keyword{models}
