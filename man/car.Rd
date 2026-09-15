\name{car}
\alias{car}
\title{Proper conditional autoregressive spatial covariance structure}
\description{
Proper conditional autoregressive spatial covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
car(x, W, rho = 0.10, fixed = FALSE)
}
\arguments{
  \item{x}{Factor or design defining q spatial levels.}
  \item{W}{Symmetric nonnegative adjacency/weights matrix with zero diagonal, aligned to the levels of x.}
  \item{rho}{Starting proper-CAR dependence parameter.}
  \item{fixed}{Logical indicating whether rho is fixed.}
}
\details{
Let \eqn{D=\mathrm{diag}(W\mathbf{1})}.  The proper CAR precision and covariance are \deqn{Q=D-\rho W,\qquad M=Q^{-1},\qquad K=M/M_{11}.} The constructor requires positive row sums, so isolated levels are not allowed.  The admissible open interval for \code{rho} is obtained from the eigenvalues of \eqn{D^{-1/2}WD^{-1/2}} and is enforced through a bounded-logit working coordinate.  This is a proper, nonsingular CAR model; an intrinsic singular CAR is not used by the current precision-based Henderson implementation.  The covariance derivative is analytic.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(car(location, W), ism(genotype))
}
}
\seealso{
\code{\link{sar}}, \code{\link{maternm}}, \code{\link{vsm}}.
}
