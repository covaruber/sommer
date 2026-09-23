\name{scores_mmes}
\alias{scores_mmes}
\title{Predict factor-analytic scores from a fitted mmes model}
\description{
Predicts latent factor scores for each level of the main random-effect term
crossed with a \code{fam()} or \code{rrcm()} covariance-shaping factor,
combining the fitted loadings (see \code{\link{loadings_mmes}}), the fitted
covariance matrix, and the BLUPs already stored in \code{object$uList}.
}
\usage{
scores_mmes(object, term = NULL, method = c("regression", "bartlett"), 
            varianceScale=TRUE, rotation=TRUE)
}
\arguments{
  \item{object}{a fitted model of class \code{"mmes"}.}
  \item{term}{character name of the random term built with a single
\code{fam()} or \code{rrcm()} covariance-shaping factor. If \code{NULL}, the
unique such term in the model is used automatically.}
  \item{method}{\code{"regression"} (default) uses the full fitted covariance
\eqn{\Sigma} (Thomson's method); \code{"bartlett"} uses only the specific
(residual) variances \eqn{\Psi} (Bartlett's classic unbiased estimator).}
\item{varianceScale}{a logical argument to indicate if loadings should be returned 
in variance scale (multiplied by sqrt(sigma2)).}
\item{rotation}{a logical value to indicate if loadings should be rotated by its
singular vectors.}
}
\details{
Let \eqn{U} be the levels-by-levels BLUP matrix in \code{object$uList[[term]]}
and \eqn{\Lambda} the loadings from \code{\link{loadings_mmes}}. The
regression (Thomson) scores are
\deqn{F = U\Sigma^{-1}\Lambda,}
and the Bartlett scores are
\deqn{F = U\Psi^{-1}\Lambda(\Lambda^{\mathsf T}\Psi^{-1}\Lambda)^{-1}.}
Because \eqn{U} already contains shrunken (BLUP) effects, the resulting
scores are similarly regularized.
}
\value{
A matrix with one row per level of the main random-effect term and one
column per latent factor.
}
\examples{
\dontrun{
mix <- mmes(BLUEs ~ trial,
            random = ~ vsm(fam(trial, 2), ism(genotype)),
            rcov = ~ units, data = dt)
scores_mmes(mix)
}
}
\seealso{
\code{\link{loadings_mmes}}, \code{\link{fam}}, \code{\link{rrcm}}, \code{\link{vsm}}, \code{\link{mmes}}.
}
